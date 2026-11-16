import Mathlib.Algebra.BigOperators.Field
import Mathlib.Algebra.BigOperators.Ring.Finset
import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Algebra.Order.BigOperators.Ring.Finset
import Mathlib.Analysis.Calculus.Deriv.MeanValue
import Mathlib.Analysis.Calculus.ParametricIntegral
import Mathlib.Analysis.SpecialFunctions.Exponential
import Mathlib.Analysis.SpecialFunctions.ImproperIntegrals
import Mathlib.Analysis.SpecialFunctions.Integrals.Basic
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Analysis.SpecialFunctions.Log.Deriv
import Mathlib.Analysis.SpecialFunctions.Sqrt
import Mathlib.Data.Finset.Card
import Mathlib.Data.Finset.Max
import Mathlib.Data.Fintype.BigOperators
import Mathlib.Data.Fintype.Perm
import Mathlib.Data.Prod.Lex
import Mathlib.Data.Real.Basic
import Mathlib.GroupTheory.Perm.Fin
import Mathlib.MeasureTheory.Constructions.Pi
import Mathlib.MeasureTheory.Group.Convolution
import Mathlib.MeasureTheory.Group.LIntegral
import Mathlib.MeasureTheory.Integral.DominatedConvergence
import Mathlib.MeasureTheory.Integral.Lebesgue.Markov
import Mathlib.MeasureTheory.Measure.Lebesgue.Basic
import Mathlib.Order.Interval.Set.OrdConnected
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Ring
import Mathlib.Topology.Order.IntermediateValue

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

/-! ### E35's follow-up: M-1 under covariate shift -/

section CovariateShift

variable {Z β : Type*} [LinearOrder β]

/-- `P(L + θT ≤ u)` for `θ` uniform on `(0, 1]`: a smoothed p-value's distribution function
given the
mass `L` strictly below its score and the mass `T` tied with it, itself included. -/
noncomputable def smoothCdf (L T u : ℝ) : ℝ :=
  if 0 < T then max 0 (min 1 ((u - L) / T)) else if L ≤ u then 1 else 0

/-- An ordering's weight, newest first: the likelihood ratio over its newest `s` positions, the
stream's, and 1 over the rest, the reference sample's (Tibshirani et al.'s Definition 1). -/
noncomputable def streamWeight (w : Z → ℝ) (s : ℕ) {n : ℕ} (x : Fin n → Z) : ℝ :=
  ∏ i : Fin n, if (i : ℕ) < s then w (x i) else 1

/-- The weights of every ordering of the bag `x`. -/
noncomputable def streamTotal (w : Z → ℝ) (s : ℕ) {n : ℕ} (x : Fin n → Z) : ℝ :=
  ∑ τ : Equiv.Perm (Fin n), streamWeight w s (x ∘ τ)

/-- **Tibshirani et al.'s Lemma 3:** the chance that unit `j` of the bag holds the newest position,
each ordering weighed by its stream weight. -/
noncomputable def newestShare (w : Z → ℝ) (s : ℕ) {n : ℕ} (x : Fin (n + 1) → Z) (j : Fin (n +
    1)) : ℝ :=
  (∑ τ ∈ univ.filter fun τ : Equiv.Perm (Fin (n + 1)) => τ 0 = j, streamWeight w s (x ∘ τ)) /
    streamTotal w s x

/-- The weighted mass strictly below the newest score. -/
noncomputable def shiftBelow (w : Z → ℝ) (a : Z → β) (s : ℕ) {n : ℕ} (x : Fin (n + 1) → Z) : ℝ :=
  ∑ j, if a (x j) < a (x 0) then newestShare w s x j else 0

/-- The weighted mass tied with the newest score, itself included. -/
noncomputable def shiftTied (w : Z → ℝ) (a : Z → β) (s : ℕ) {n : ℕ} (x : Fin (n + 1) → Z) : ℝ :=
  ∑ j, if a (x j) = a (x 0) then newestShare w s x j else 0

/-- Newest first: the product over the newest `k` arrivals, the stream's, of each weighted smoothed
p-value's distribution function, at `us 0` for the newest. At each arrival the bag is every unit so
far, its stream positions the newest; the oldest `n − k` units are the reference sample. -/
noncomputable def shiftSeqCdf (w : Z → ℝ) (a : Z → β) :
    (k : ℕ) → {n : ℕ} → (Fin n → Z) → (Fin k → ℝ) → ℝ
  | 0, _, _, _ => 1
  | _ + 1, 0, _, _ => 0
  | k + 1, _ + 1, x, us => smoothCdf (shiftBelow w a (k + 1) x) (shiftTied w a (k + 1) x) (us 0) *
      shiftSeqCdf w a k (fun i => x i.succ) (fun i => us i.succ)

/-- The weighted tie classes cover `[0, min(u, total)]`: with masses `P ≥ 0`, each unit's smoothed
distribution function at `u ≥ 0`, weighed by its mass, sums to `min(u, ∑ P)`. -/
theorem sum_smoothCdf {ι : Type*} (P : ι → ℝ) (hP : ∀ i, 0 ≤ P i) (a : ι → β)
    {u : ℝ} (hu : 0 ≤ u) : ∀ C : Finset ι,
      ∑ j ∈ C, P j * smoothCdf (∑ i ∈ C, if a i < a j then P i else 0)
        (∑ i ∈ C, if a i = a j then P i else 0) u = min u (∑ i ∈ C, P i) := by
  classical
  intro C
  induction C using Finset.strongInduction with
  | H C ih =>
    rcases C.eq_empty_or_nonempty with hC | hC
    · subst hC; simp [hu]
    obtain ⟨t, ht, htmax⟩ := C.exists_max_image a hC
    set M := C.filter fun j => a j = a t with hM
    set C' := C.filter fun j => a j ≠ a t with hC'
    have hsplit : ∀ f : ι → ℝ, ∑ j ∈ C, f j = ∑ j ∈ C', f j + ∑ j ∈ M, f j := fun f => by
      rw [hC', hM, add_comm, sum_filter_add_sum_filter_not]
    have hsub : C' ⊂ C := by
      refine (filter_ssubset).2 ⟨t, ht, by simp⟩
    have ih' := ih C' hsub
    -- below the top class, the masses read in C are those in C'
    have hlow : ∀ j ∈ C', (∑ i ∈ C, if a i < a j then P i else 0) =
        ∑ i ∈ C', if a i < a j then P i else 0 := fun j hj => by
      rw [hsplit]
      have : ∑ i ∈ M, (if a i < a j then P i else 0) = 0 := sum_eq_zero fun i hi => by
        have hi' := (mem_filter.1 hi).2
        have hj' := mem_filter.1 hj
        have : ¬ a i < a j := by rw [hi']; exact not_lt.2 (htmax j hj'.1)
        simp [this]
      rw [this, add_zero]
    have htie : ∀ j ∈ C', (∑ i ∈ C, if a i = a j then P i else 0) =
        ∑ i ∈ C', if a i = a j then P i else 0 := fun j hj => by
      rw [hsplit]
      have : ∑ i ∈ M, (if a i = a j then P i else 0) = 0 := sum_eq_zero fun i hi => by
        have hi' := (mem_filter.1 hi).2
        have hj' := (mem_filter.1 hj).2
        have : ¬ a i = a j := by rw [hi']; exact fun h => hj' h.symm
        simp [this]
      rw [this, add_zero]
    -- the top class reads everything else below it, and itself tied
    set S' := ∑ i ∈ C', P i
    set TM := ∑ i ∈ M, P i
    have htopL : ∀ j ∈ M, (∑ i ∈ C, if a i < a j then P i else 0) = S' := fun j hj => by
      have hj' := (mem_filter.1 hj).2
      rw [hsplit]
      have h1 : ∑ i ∈ C', (if a i < a j then P i else 0) = S' := sum_congr rfl fun i hi => by
        have hi' := mem_filter.1 hi
        have : a i < a j := by rw [hj']; exact lt_of_le_of_ne (htmax i hi'.1) hi'.2
        simp [this]
      have h2 : ∑ i ∈ M, (if a i < a j then P i else 0) = 0 := sum_eq_zero fun i hi => by
        have hi' := (mem_filter.1 hi).2
        simp [hi', hj']
      rw [h1, h2, add_zero]
    have htopT : ∀ j ∈ M, (∑ i ∈ C, if a i = a j then P i else 0) = TM := fun j hj => by
      have hj' := (mem_filter.1 hj).2
      rw [hsplit]
      have h1 : ∑ i ∈ C', (if a i = a j then P i else 0) = 0 := sum_eq_zero fun i hi => by
        have hi' := (mem_filter.1 hi).2
        have : ¬ a i = a j := by rw [hj']; exact hi'
        simp [this]
      have h2 : ∑ i ∈ M, (if a i = a j then P i else 0) = TM := sum_congr rfl fun i hi => by
        have hi' := (mem_filter.1 hi).2
        simp [hi', hj']
      rw [h1, h2, zero_add]
    have hS'0 : 0 ≤ S' := sum_nonneg fun i _ => hP i
    have hTM0 : 0 ≤ TM := sum_nonneg fun i _ => hP i
    rw [hsplit, hsplit P]
    rw [sum_congr rfl fun j hj => by rw [hlow j hj, htie j hj], ih']
    rw [sum_congr rfl fun j hj => by rw [htopL j hj, htopT j hj], ← sum_mul]
    change min u S' + TM * smoothCdf S' TM u = min u (S' + TM)
    unfold smoothCdf
    rcases eq_or_lt_of_le hTM0 with h0 | hpos
    · rw [← h0]; simp
    · rw [if_pos hpos]
      have key : TM * max 0 (min 1 ((u - S') / TM)) = max 0 (min TM (u - S')) := by
        rw [mul_max_of_nonneg _ _ hTM0, mul_min_of_nonneg _ _ hTM0, mul_zero, mul_one,
          mul_div_cancel₀ _ hpos.ne']
      rw [key]
      simp only [min_def, max_def]
      split_ifs <;> linarith

/-- The bag without unit `p`, in the order `decomposeFin` gives it. -/
def tailAt {n : ℕ} (x : Fin (n + 1) → Z) (p : Fin (n + 1)) : Fin n → Z :=
  fun i => x (Equiv.swap 0 p i.succ)

theorem comp_decomposeFin_zero {n : ℕ} (x : Fin (n + 1) → Z) (p : Fin (n + 1))
    (e : Equiv.Perm (Fin n)) : (x ∘ Equiv.Perm.decomposeFin.symm (p, e)) 0 = x p := by
  simp [Equiv.Perm.decomposeFin_symm_apply_zero]

theorem comp_decomposeFin_succ {n : ℕ} (x : Fin (n + 1) → Z) (p : Fin (n + 1))
    (e : Equiv.Perm (Fin n)) :
    (fun i : Fin n => (x ∘ Equiv.Perm.decomposeFin.symm (p, e)) i.succ) = tailAt x p ∘ e := by
  funext i
  simp [tailAt, Equiv.Perm.decomposeFin_symm_apply_succ]

theorem streamWeight_succ (w : Z → ℝ) (s : ℕ) {n : ℕ} (x : Fin (n + 1) → Z) :
    streamWeight w (s + 1) x = w (x 0) * streamWeight w s (fun i => x i.succ) := by
  unfold streamWeight
  rw [Fin.prod_univ_succ]
  simp [Fin.val_succ]

theorem streamWeight_nonneg {w : Z → ℝ} (hw : ∀ z, 0 ≤ w z) (s : ℕ) {n : ℕ} (x : Fin n → Z) :
    0 ≤ streamWeight w s x :=
  prod_nonneg fun i _ => by split_ifs; exacts [hw _, zero_le_one]

theorem streamTotal_nonneg {w : Z → ℝ} (hw : ∀ z, 0 ≤ w z) (s : ℕ) {n : ℕ} (x : Fin n → Z) :
    0 ≤ streamTotal w s x :=
  sum_nonneg fun _ _ => streamWeight_nonneg hw s _

/-- Placing a unit at the newest position: the bag's weight is the newest's ratio times the
rest's. -/
theorem streamTotal_succ (w : Z → ℝ) (s : ℕ) {n : ℕ} (x : Fin (n + 1) → Z) :
    streamTotal w (s + 1) x = ∑ p, w (x p) * streamTotal w s (tailAt x p) := by
  unfold streamTotal
  rw [← Equiv.sum_comp Equiv.Perm.decomposeFin.symm, Fintype.sum_prod_type]
  refine sum_congr rfl fun p _ => ?_
  rw [mul_sum]
  refine sum_congr rfl fun e _ => ?_
  rw [streamWeight_succ, comp_decomposeFin_zero, comp_decomposeFin_succ]

/-- The orderings putting unit `p` newest weigh its ratio times the rest's. -/
theorem sum_newest_eq (w : Z → ℝ) (s : ℕ) {n : ℕ} (x : Fin (n + 1) → Z) (p : Fin (n + 1)) :
    ∑ τ ∈ univ.filter fun τ : Equiv.Perm (Fin (n + 1)) => τ 0 = p, streamWeight w (s + 1) (x ∘ τ) =
      w (x p) * streamTotal w s (tailAt x p) := by
  rw [sum_filter, ← Equiv.sum_comp Equiv.Perm.decomposeFin.symm, Fintype.sum_prod_type]
  simp only [Equiv.Perm.decomposeFin_symm_apply_zero]
  rw [sum_eq_single p]
  · simp only [if_true]
    unfold streamTotal
    rw [mul_sum]
    refine sum_congr rfl fun e _ => ?_
    rw [streamWeight_succ, comp_decomposeFin_zero, comp_decomposeFin_succ]
  · intro q _ hq
    exact sum_eq_zero fun e _ => if_neg hq
  · intro h; exact absurd (mem_univ p) h

theorem streamTotal_comp (w : Z → ℝ) (s : ℕ) {n : ℕ} (x : Fin n → Z) (σ : Equiv.Perm (Fin n)) :
    streamTotal w s (x ∘ σ) = streamTotal w s x := by
  unfold streamTotal
  exact Equiv.sum_comp (Equiv.mulLeft σ) (fun τ => streamWeight w s (x ∘ τ))

theorem newestShare_comp (w : Z → ℝ) (s : ℕ) {n : ℕ} (x : Fin (n + 1) → Z)
    (σ : Equiv.Perm (Fin (n + 1))) (j : Fin (n + 1)) :
    newestShare w s (x ∘ σ) j = newestShare w s x (σ j) := by
  unfold newestShare
  rw [streamTotal_comp]
  congr 1
  rw [sum_filter, sum_filter]
  rw [← Equiv.sum_comp (Equiv.mulLeft σ⁻¹)]
  refine sum_congr rfl fun τ _ => ?_
  simp only [Equiv.coe_mulLeft, Equiv.Perm.coe_mul, Function.comp_apply]
  by_cases h : τ 0 = σ j
  · have h' : σ⁻¹ (τ 0) = j := by rw [h]; simp
    rw [if_pos h', if_pos h]
    congr 1
    funext i
    simp
  · have h' : ¬ σ⁻¹ (τ 0) = j := fun hh => h (by rw [← hh]; simp)
    rw [if_neg h', if_neg h]

theorem shiftBelow_comp (w : Z → ℝ) (a : Z → β) (s : ℕ) {n : ℕ} (x : Fin (n + 1) → Z)
    (σ : Equiv.Perm (Fin (n + 1))) :
    shiftBelow w a s (x ∘ σ) = ∑ j, if a (x j) < a (x (σ 0)) then newestShare w s x j else 0 := by
  unfold shiftBelow
  simp only [newestShare_comp, Function.comp_apply]
  exact Equiv.sum_comp σ (fun j => if a (x j) < a (x (σ 0)) then newestShare w s x j else 0)

theorem shiftTied_comp (w : Z → ℝ) (a : Z → β) (s : ℕ) {n : ℕ} (x : Fin (n + 1) → Z)
    (σ : Equiv.Perm (Fin (n + 1))) :
    shiftTied w a s (x ∘ σ) = ∑ j, if a (x j) = a (x (σ 0)) then newestShare w s x j else 0 := by
  unfold shiftTied
  simp only [newestShare_comp, Function.comp_apply]
  exact Equiv.sum_comp σ (fun j => if a (x j) = a (x (σ 0)) then newestShare w s x j else 0)

/-- **M-1 under covariate shift, the finite core.** Over every ordering of a fixed bag of `n`
units, each
weighed by the likelihood ratio over its newest `k` positions, the stream's (Tibshirani et al.'s
Definition 1, the oldest `n − k` units the reference sample), the stream's weighted smoothed
p-values
have joint distribution function `∏ uₜ`: independent and uniform. Each p-value ranks its unit among
every unit so far by Lemma 3's chance of holding the newest position, its ties broken by a fresh
uniform draw. -/
theorem shift_uniform (w : Z → ℝ) (hw : ∀ z, 0 ≤ w z) (a : Z → β) :
    ∀ (k : ℕ) {n : ℕ} (x : Fin n → Z) (us : Fin k → ℝ), k ≤ n → (∀ i, 0 ≤ us i ∧ us i ≤ 1) →
      ∑ σ : Equiv.Perm (Fin n), streamWeight w k (x ∘ σ) * shiftSeqCdf w a k (x ∘ σ) us =
        streamTotal w k x * ∏ i, us i
  | 0, n, x, us, _, _ => by simp [shiftSeqCdf, streamTotal]
  | k + 1, 0, _, _, hk, _ => absurd hk (by omega)
  | k + 1, n + 1, x, us, hk, hus => by
    have ih := fun y : Fin n → Z =>
      shift_uniform w hw a k y (fun i => us i.succ) (by omega) fun i => hus i.succ
    set G : Fin (n + 1) → ℝ := fun p =>
      smoothCdf (∑ j, if a (x j) < a (x p) then newestShare w (k + 1) x j else 0)
        (∑ j, if a (x j) = a (x p) then newestShare w (k + 1) x j else 0) (us 0) with hG
    have hstep : ∀ (p : Fin (n + 1)) (e : Equiv.Perm (Fin n)),
        streamWeight w (k + 1) (x ∘ Equiv.Perm.decomposeFin.symm (p, e)) *
            shiftSeqCdf w a (k + 1) (x ∘ Equiv.Perm.decomposeFin.symm (p, e)) us =
          (w (x p) * G p) * (streamWeight w k (tailAt x p ∘ e) *
            shiftSeqCdf w a k (tailAt x p ∘ e) (fun i => us i.succ)) := by
      intro p e
      have hz := comp_decomposeFin_zero x p e
      have hs := comp_decomposeFin_succ x p e
      simp only [shiftSeqCdf]
      rw [streamWeight_succ, shiftBelow_comp, shiftTied_comp,
          Equiv.Perm.decomposeFin_symm_apply_zero, hz,
        hs, hG]
      ring
    rw [← Equiv.sum_comp Equiv.Perm.decomposeFin.symm, Fintype.sum_prod_type]
    simp only [hstep, ← mul_sum, ih]
    have hT := streamTotal_succ w k x
    have hnn : ∀ p, 0 ≤ w (x p) * streamTotal w k (tailAt x p) := fun p =>
      mul_nonneg (hw _) (streamTotal_nonneg hw k _)
    rw [Fin.prod_univ_succ]
    rcases eq_or_lt_of_le (streamTotal_nonneg hw (k + 1) x) with h0 | hpos
    · have hz : ∀ p, w (x p) * streamTotal w k (tailAt x p) = 0 := fun p =>
        (sum_eq_zero_iff_of_nonneg fun p _ => hnn p).1 (hT ▸ h0.symm) p (mem_univ p)
      rw [← h0, zero_mul]
      exact sum_eq_zero fun p _ => by
        rw [show w (x p) * G p * (streamTotal w k (tailAt x p) * ∏ i : Fin k, us i.succ) =
          (w (x p) * streamTotal w k (tailAt x p)) * (G p * ∏ i : Fin k, us i.succ) by ring, hz
              p, zero_mul]
    · have hshare : ∀ p, newestShare w (k + 1) x p =
          w (x p) * streamTotal w k (tailAt x p) / streamTotal w (k + 1) x := fun p => by
        unfold newestShare; rw [sum_newest_eq]
      have hsum1 : ∑ p, newestShare w (k + 1) x p = 1 := by
        simp only [hshare, ← sum_div, ← hT]
        exact div_self hpos.ne'
      have hcov := sum_smoothCdf (fun p => newestShare w (k + 1) x p)
        (fun p => by rw [hshare]; exact div_nonneg (hnn p) hpos.le) (fun p => a (x p)) (hus 0).1
            univ
      simp only [hsum1, min_eq_left (hus 0).2] at hcov
      have : ∑ p, w (x p) * G p * (streamTotal w k (tailAt x p) * ∏ i : Fin k, us i.succ) =
          streamTotal w (k + 1) x * (∑ p, newestShare w (k + 1) x p * G p) * ∏ i : Fin k, us
              i.succ := by
        rw [mul_sum, sum_mul]
        refine sum_congr rfl fun p _ => ?_
        rw [hshare]
        field_simp
      rw [this, hcov]
      ring

/-- **M-1 under covariate shift, under a weighted exchangeable law** (Tibshirani et al.'s
Definition 1): a law on bags of `n` units with probability `∏ (ratios over the newest k) · g`, `g`
symmetric. Then the stream's weighted smoothed p-values have joint distribution function `∏ uₜ`. -/
theorem shift_law [Fintype Z] (w : Z → ℝ) (hw : ∀ z, 0 ≤ w z) (a : Z → β) {n k : ℕ} (hk : k ≤ n)
    (g : (Fin n → Z) → ℝ) (hg : ∀ σ : Equiv.Perm (Fin n), ∀ z, g (z ∘ σ) = g z)
    (h1 : ∑ z, streamWeight w k z * g z = 1) (us : Fin k → ℝ) (hus : ∀ i, 0 ≤ us i ∧ us i ≤ 1) :
    ∑ z, streamWeight w k z * g z * shiftSeqCdf w a k z us = ∏ i, us i := by
  have huniv : permsOn (univ : Finset (Fin n)) = univ :=
    eq_univ_of_forall fun σ => mem_permsOn_univ σ
  have hA := sum_law_average (C := univ) g (fun σ _ z => hg σ z)
    (fun z => streamWeight w k z * shiftSeqCdf w a k z us)
  have hB := sum_law_average (C := univ) g (fun σ _ z => hg σ z) (fun z => streamWeight w k z)
  rw [huniv] at hA hB
  simp only [shift_uniform w hw a k _ us hk hus] at hA
  have hB' : ∑ z, g z * (streamTotal w k z / (univ : Finset (Equiv.Perm (Fin n))).card) = 1 := by
    rw [← h1]
    unfold streamTotal
    rw [show ∑ z, streamWeight w k z * g z = ∑ z, g z * streamWeight w k z from
      sum_congr rfl fun z _ => mul_comm _ _, hB]
  calc ∑ z, streamWeight w k z * g z * shiftSeqCdf w a k z us
      = ∑ z, g z * (streamWeight w k z * shiftSeqCdf w a k z us) :=
        sum_congr rfl fun z _ => by ring
    _ = ∑ z, g z * ((streamTotal w k z * ∏ i, us i) /
          (univ : Finset (Equiv.Perm (Fin n))).card) := hA
    _ = (∑ z, g z * (streamTotal w k z / (univ : Finset (Equiv.Perm (Fin n))).card)) *
          ∏ i, us i := by
        rw [sum_mul]; exact sum_congr rfl fun z _ => by ring
    _ = ∏ i, us i := by rw [hB', one_mul]

/-- The elementary symmetric polynomial `e_r` of the values `v` over `C`. -/
noncomputable def esymmOn {ι : Type*} (C : Finset ι) (r : ℕ) (v : ι → ℝ) : ℝ :=
  ∑ S ∈ C.powersetCard r, ∏ i ∈ S, v i

theorem esymmOn_insert {ι : Type*} [DecidableEq ι] {a : ι} {D : Finset ι} (ha : a ∉ D) (r : ℕ)
    (v : ι → ℝ) : esymmOn (insert a D) (r + 1) v = esymmOn D (r + 1) v + v a * esymmOn D r v := by
  unfold esymmOn
  rw [powersetCard_succ_insert ha, sum_union]
  · congr 1
    rw [sum_image, mul_sum]
    · refine sum_congr rfl fun S hS => ?_
      rw [prod_insert fun h => ha ((mem_powersetCard.1 hS).1 h)]
    · intro S hS T hT hST
      have hS' := (mem_powersetCard.1 hS).1
      have hT' := (mem_powersetCard.1 hT).1
      rw [← erase_insert (fun h => ha (hS' h)), ← erase_insert (fun h => ha (hT' h))]
      rw [hST]
  · rw [disjoint_left]
    intro S hS hS2
    obtain ⟨T, -, rfl⟩ := mem_image.1 hS2
    exact ha ((mem_powersetCard.1 hS).1 (mem_insert_self a T))

/-- Counting each `(r + 1)`-set once per member: `∑ⱼ vⱼ e_r(v₋ⱼ) = (r + 1) e_(r+1)(v)`. -/
theorem sum_mul_esymmOn_erase {ι : Type*} [DecidableEq ι] (v : ι → ℝ) :
    ∀ (C : Finset ι) (r : ℕ), ∑ p ∈ C, v p * esymmOn (C.erase p) r v = (r + 1) * esymmOn C (r +
        1) v := by
  intro C
  induction C using Finset.induction_on with
  | empty =>
    intro r
    simp only [sum_empty, esymmOn]
    rw [powersetCard_eq_empty.2 (by simp), sum_empty, mul_zero]
  | insert a B ha ih =>
    intro r
    rw [sum_insert ha, erase_insert ha]
    have herase : ∀ p ∈ B, (insert a B).erase p = insert a (B.erase p) := fun p hp =>
      erase_insert_of_ne (fun h => ha (h ▸ hp))
    rw [sum_congr rfl fun p hp => by rw [herase p hp]]
    cases r with
    | zero =>
      simp only [esymmOn, powersetCard_zero, sum_singleton, prod_empty, mul_one, Nat.cast_zero,
        zero_add, one_mul, powersetCard_one, sum_map, Function.Embedding.coeFn_mk, prod_singleton]
      rw [sum_insert ha]
    | succ r =>
      have hins : ∀ p ∈ B, esymmOn (insert a (B.erase p)) (r + 1) v =
          esymmOn (B.erase p) (r + 1) v + v a * esymmOn (B.erase p) r v := fun p hp =>
        esymmOn_insert (fun h => ha (mem_of_mem_erase h)) r v
      rw [sum_congr rfl fun p hp => by rw [hins p hp]]
      simp only [mul_add, sum_add_distrib]
      rw [ih (r + 1), esymmOn_insert ha (r + 1) v]
      have h2 : ∑ p ∈ B, v p * (v a * esymmOn (B.erase p) r v) =
          v a * ∑ p ∈ B, v p * esymmOn (B.erase p) r v := by
        rw [mul_sum]; exact sum_congr rfl fun p _ => by ring
      rw [h2, ih r]
      push_cast
      ring

theorem esymmOn_map {ι κ : Type*} (C : Finset ι) (f : ι ↪ κ) (r : ℕ) (v : κ → ℝ) :
    esymmOn (C.map f) r v = esymmOn C r (v ∘ f) := by
  unfold esymmOn
  rw [powersetCard_map, sum_map]
  refine sum_congr rfl fun S _ => ?_
  simp [prod_map]

/-- The bag without unit `p` is the bag's other units. -/
theorem map_tailAt (n : ℕ) (p : Fin (n + 1)) :
    (univ : Finset (Fin n)).map ⟨fun i => Equiv.swap 0 p i.succ, fun i j h => by
      simpa using h⟩ = univ.erase p := by
  apply eq_of_subset_of_card_le
  · intro q hq
    obtain ⟨i, -, rfl⟩ := mem_map.1 hq
    refine mem_erase.2 ⟨fun h => ?_, mem_univ _⟩
    simp only [Function.Embedding.coeFn_mk] at h
    have : (i.succ : Fin (n + 1)) = 0 := by
      have h' := congrArg (Equiv.swap 0 p) h
      simp at h'
    exact Fin.succ_ne_zero i this
  · rw [card_map, card_erase_of_mem (mem_univ p), card_univ, card_univ, Fintype.card_fin,
      Fintype.card_fin]
    omega

/-- **The orderings' total weight is `s! (n − s)! e_s(w)`**: each set of `s` units fills the stream
positions in `s!` orders and the reference positions in `(n − s)!`. -/
theorem streamTotal_esymm (w : Z → ℝ) :
    ∀ (s : ℕ) {n : ℕ} (x : Fin n → Z), s ≤ n →
      streamTotal w s x = s.factorial * (n - s).factorial * esymmOn univ s (w ∘ x)
  | 0, n, x, _ => by
    simp [streamTotal, streamWeight, esymmOn, Fintype.card_perm]
  | s + 1, 0, _, hs => absurd hs (by omega)
  | s + 1, n + 1, x, hs => by
    rw [streamTotal_succ]
    have hmap : ∀ p : Fin (n + 1), esymmOn univ s (w ∘ tailAt x p) =
        esymmOn ((univ : Finset (Fin (n + 1))).erase p) s (w ∘ x) := fun p => by
      rw [← map_tailAt n p, esymmOn_map]
      rfl
    simp only [streamTotal_esymm w s _ (by omega : s ≤ n), hmap]
    have hsum := sum_mul_esymmOn_erase (w ∘ x) univ s
    simp only [Function.comp_apply] at hsum
    rw [show ∑ p : Fin (n + 1), w (x p) * ((s.factorial : ℝ) * ((n - s).factorial : ℝ) *
        esymmOn (univ.erase p) s (w ∘ x)) = (s.factorial : ℝ) * ((n - s).factorial : ℝ) *
        ∑ p : Fin (n + 1), w (x p) * esymmOn (univ.erase p) s (w ∘ x) by
      rw [mul_sum]; exact sum_congr rfl fun p _ => by ring]
    rw [hsum, show n + 1 - (s + 1) = n - s by omega, Nat.factorial_succ]
    push_cast
    ring

/-- **The share, as `CovariateShiftPValues` computes it**: Lemma 3's chance that unit `j` holds the
newest of `k + 1` stream positions is `πⱼ/(k + 1)`, with `πⱼ = wⱼ e_k(w₋ⱼ)/e_(k+1)(w)` the
conditional
Poisson inclusion probability. -/
theorem newestShare_eq_esymm (w : Z → ℝ) {k n : ℕ} (hk : k ≤ n) (x : Fin (n + 1) → Z) (j : Fin
    (n + 1)) :
    newestShare w (k + 1) x j = w (x j) * esymmOn (univ.erase j) k (w ∘ x) /
      ((k + 1) * esymmOn univ (k + 1) (w ∘ x)) := by
  unfold newestShare
  rw [sum_newest_eq, streamTotal_esymm w k _ hk, streamTotal_esymm w (k + 1) _ (by omega)]
  have hmap : esymmOn univ k (w ∘ tailAt x j) = esymmOn ((univ : Finset (Fin (n + 1))).erase j)
      k (w ∘ x) := by
    rw [← map_tailAt n j, esymmOn_map]
    rfl
  rw [hmap, show n + 1 - (k + 1) = n - k by omega, Nat.factorial_succ]
  have hf : (0 : ℝ) < k.factorial := by exact_mod_cast k.factorial_pos
  have hg : (0 : ℝ) < (n - k).factorial := by exact_mod_cast (n - k).factorial_pos
  push_cast
  rcases eq_or_ne (esymmOn univ (k + 1) (w ∘ x)) 0 with h0 | h0
  · simp [h0]
  · field_simp

end CovariateShift

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

/-! ### E35's follow-up: CUSUM's and Shiryaev–Roberts's guarantees -/

section ChangeAlarms

open MeasureTheory
open scoped ENNReal

variable (μ : Measure ℝ) (f : (n : ℕ) → (Fin n → ℝ) → ℝ → ℝ≥0∞) (c : ℝ≥0∞)

/-- One Shiryaev–Roberts step on the bet's value `b` from `(Q, Σ, A)`: the statistic
`R = b (Q + 1)`; at `c ≤ R` an alarm, which adds `R` to `Σ`, counts in `A`, and restarts `Q` at
0. -/
noncomputable def srStep (b : ℝ≥0∞) (s : ℝ≥0∞ × ℝ≥0∞ × ℕ) : ℝ≥0∞ × ℝ≥0∞ × ℕ :=
  if c ≤ b * (s.1 + 1) then (0, s.2.1 + b * (s.1 + 1), s.2.2 + 1) else (b * (s.1 + 1), s.2.1, s.2.2)

/-- Shiryaev–Roberts with restarts (`MartingaleAlarm`'s): `Rₙ = fₙ (Rₙ₋₁ + 1)` from `R₀ = 0`, the
`n`th bet chosen from the earlier p-values and read at the newest. -/
noncomputable def srRun : (n : ℕ) → (Fin n → ℝ) → ℝ≥0∞ × ℝ≥0∞ × ℕ
  | 0, _ => (0, 0, 0)
  | n + 1, h => srStep c (f n (Fin.init h) (h (Fin.last n))) (srRun n (Fin.init h))

/-- One CUSUM step from `(γ, A)`: Page's `γ = b max(γ, 1)`; at `c ≤ γ` an alarm and a restart. -/
noncomputable def cusumStep (b : ℝ≥0∞) (s : ℝ≥0∞ × ℕ) : ℝ≥0∞ × ℕ :=
  if c ≤ b * max s.1 1 then (0, s.2 + 1) else (b * max s.1 1, s.2)

/-- CUSUM with restarts (`MartingaleAlarm`'s): `γₙ = fₙ max(γₙ₋₁, 1)` from `γ₀ = 0`. -/
noncomputable def cusumRun : (n : ℕ) → (Fin n → ℝ) → ℝ≥0∞ × ℕ
  | 0, _ => (0, 0)
  | n + 1, h => cusumStep c (f n (Fin.init h) (h (Fin.last n))) (cusumRun n (Fin.init h))

theorem srStep_sum (b : ℝ≥0∞) (s : ℝ≥0∞ × ℝ≥0∞ × ℕ) :
    (srStep c b s).1 + (srStep c b s).2.1 = b * (s.1 + 1) + s.2.1 := by
  unfold srStep
  split_ifs <;> simp [add_comm]

theorem srStep_count (b : ℝ≥0∞) (s : ℝ≥0∞ × ℝ≥0∞ × ℕ) (hs : c * s.2.2 ≤ s.2.1) :
    c * (srStep c b s).2.2 ≤ (srStep c b s).2.1 := by
  unfold srStep
  split_ifs with h
  · dsimp only
    have e : c * (((s.2.2 + 1 : ℕ)) : ℝ≥0∞) = c * s.2.2 + c := by push_cast; rw [mul_add, mul_one]
    rw [e]
    exact add_le_add hs h
  · exact hs

theorem srRun_count : ∀ (n : ℕ) (h : Fin n → ℝ), c * (srRun f c n h).2.2 ≤ (srRun f c n h).2.1
  | 0, _ => by simp [srRun]
  | n + 1, h => srStep_count c _ _ (srRun_count n (Fin.init h))

theorem measurable_srStep {α : Type*} [MeasurableSpace α] {b : α → ℝ≥0∞}
    {s : α → ℝ≥0∞ × ℝ≥0∞ × ℕ} (hb : Measurable b) (hs : Measurable s) :
    Measurable fun a => srStep c (b a) (s a) := by
  unfold srStep
  have hR : Measurable fun a => b a * ((s a).1 + 1) := hb.mul (hs.fst.add measurable_const)
  refine Measurable.ite (measurableSet_le measurable_const hR) ?_ ?_
  · exact measurable_const.prodMk ((hs.snd.fst.add hR).prodMk
      ((measurable_from_top (f := fun k : ℕ => k + 1)).comp hs.snd.snd))
  · exact hR.prodMk (hs.snd.fst.prodMk hs.snd.snd)

variable {f} in
theorem measurable_srRun (hf : ∀ n, Measurable (Function.uncurry (f n))) :
    ∀ n, Measurable (srRun f c n)
  | 0 => measurable_const
  | n + 1 => by
    have hb := (hf n).comp ((measurable_init n).prodMk (measurable_pi_apply (Fin.last n)))
    exact measurable_srStep c hb ((measurable_srRun hf n).comp (measurable_init n))

/-- Tonelli at the newest p-value, for any measurable function of the history. -/
theorem lintegral_pi_succ [IsProbabilityMeasure μ] (n : ℕ) {g : (Fin (n + 1) → ℝ) → ℝ≥0∞}
    (hg : Measurable g) :
    ∫⁻ h, g h ∂(Measure.pi fun _ : Fin (n + 1) => μ) =
      ∫⁻ h, ∫⁻ x, g (Fin.snoc h x) ∂μ ∂(Measure.pi fun _ : Fin n => μ) := by
  set e := MeasurableEquiv.piFinSuccAbove (fun _ : Fin (n + 1) => ℝ) (Fin.last n)
  have hmp : MeasurePreserving e (Measure.pi fun _ : Fin (n + 1) => μ)
      (μ.prod (Measure.pi fun _ : Fin n => μ)) :=
    measurePreserving_piFinSuccAbove (fun _ : Fin (n + 1) => μ) (Fin.last n)
  have hG : Measurable fun z : ℝ × (Fin n → ℝ) => g (e.symm z) := hg.comp e.symm.measurable
  have hsnoc : ∀ (h : Fin n → ℝ) (x : ℝ), e.symm (x, h) = Fin.snoc h x := fun h x => by
    ext i
    simp [e, MeasurableEquiv.piFinSuccAbove_symm_apply]
  calc ∫⁻ h, g h ∂(Measure.pi fun _ : Fin (n + 1) => μ)
      = ∫⁻ h, g (e.symm (e h)) ∂(Measure.pi fun _ : Fin (n + 1) => μ) := by
        simp only [MeasurableEquiv.symm_apply_apply]
    _ = ∫⁻ z, g (e.symm z) ∂(μ.prod (Measure.pi fun _ : Fin n => μ)) :=
        hmp.lintegral_comp_emb e.measurableEmbedding fun z => g (e.symm z)
    _ = ∫⁻ h, ∫⁻ x, g (e.symm (x, h)) ∂μ ∂(Measure.pi fun _ : Fin n => μ) :=
        lintegral_prod_symm _ hG.aemeasurable
    _ = ∫⁻ h, ∫⁻ x, g (Fin.snoc h x) ∂μ ∂(Measure.pi fun _ : Fin n => μ) := by
        simp only [hsnoc]

variable {f} in
/-- **Shiryaev–Roberts's statistic and its alarms grow by at most 1 a step.** With the p-values
independent with law `μ` and each bet chosen from the past and integrating to at most 1,
`E[Qₙ + Σₙ] ≤ n`: each step adds `fₙ (Q + 1)`, which integrates to at most `Q + 1`. -/
theorem lintegral_srRun_le [IsProbabilityMeasure μ] (hf : ∀ n, Measurable (Function.uncurry (f n)))
    (h1 : ∀ n h, ∫⁻ x, f n h x ∂μ ≤ 1) :
    ∀ n, ∫⁻ h, (srRun f c n h).1 + (srRun f c n h).2.1 ∂(Measure.pi fun _ : Fin n => μ) ≤ n
  | 0 => by simp [srRun]
  | n + 1 => by
    have hm := measurable_srRun c hf
    rw [lintegral_pi_succ μ n (g := fun h => (srRun f c (n + 1) h).1 + (srRun f c (n + 1) h).2.1)
      ((hm (n + 1)).fst.add (hm (n + 1)).snd.fst)]
    have step : ∀ h : Fin n → ℝ, ∫⁻ x, (srRun f c (n + 1) (Fin.snoc h x)).1 +
        (srRun f c (n + 1) (Fin.snoc h x)).2.1 ∂μ ≤
        (srRun f c n h).1 + (srRun f c n h).2.1 + 1 := by
      intro h
      have hfx : Measurable fun x => f n h x := (hf n).comp (measurable_const.prodMk measurable_id)
      simp only [srRun, Fin.init_snoc, Fin.snoc_last, srStep_sum]
      rw [lintegral_add_right _ measurable_const, lintegral_mul_const _ hfx, lintegral_const,
        measure_univ, mul_one]
      calc (∫⁻ x, f n h x ∂μ) * ((srRun f c n h).1 + 1) + (srRun f c n h).2.1
          ≤ 1 * ((srRun f c n h).1 + 1) + (srRun f c n h).2.1 := by gcongr; exact h1 n h
        _ = (srRun f c n h).1 + (srRun f c n h).2.1 + 1 := by rw [one_mul]; ring
    calc ∫⁻ h, ∫⁻ x, (srRun f c (n + 1) (Fin.snoc h x)).1 +
          (srRun f c (n + 1) (Fin.snoc h x)).2.1 ∂μ ∂(Measure.pi fun _ : Fin n => μ)
        ≤ ∫⁻ h, (srRun f c n h).1 + (srRun f c n h).2.1 + 1 ∂(Measure.pi fun _ : Fin n => μ) :=
          lintegral_mono step
      _ = ∫⁻ h, (srRun f c n h).1 + (srRun f c n h).2.1 ∂(Measure.pi fun _ : Fin n => μ) + 1 := by
          rw [lintegral_add_right _ measurable_const, lintegral_const, measure_univ, mul_one]
      _ ≤ n + 1 := by gcongr; exact lintegral_srRun_le hf h1 n
      _ = ((n + 1 : ℕ) : ℝ≥0∞) := by push_cast; ring

variable {f} in
/-- **Shiryaev–Roberts's false alarms** (Vovk 2021, Corollary 4.6, in expectation at every
horizon): under the
null,
the expected number of alarms among the first `N` arrivals is at most `N/c`. -/
theorem sr_alarms_le [IsProbabilityMeasure μ] (hf : ∀ n, Measurable (Function.uncurry (f n)))
    (h1 : ∀ n h, ∫⁻ x, f n h x ∂μ ≤ 1) (N : ℕ) :
    c * ∫⁻ h, ((srRun f c N h).2.2 : ℝ≥0∞) ∂(Measure.pi fun _ : Fin N => μ) ≤ N := by
  have hA : Measurable fun h => ((srRun f c N h).2.2 : ℝ≥0∞) :=
    (measurable_from_top (f := fun k : ℕ => (k : ℝ≥0∞))).comp (measurable_srRun c hf N).snd.snd
  rw [← lintegral_const_mul _ hA]
  calc ∫⁻ h, c * ((srRun f c N h).2.2 : ℝ≥0∞) ∂(Measure.pi fun _ : Fin N => μ)
      ≤ ∫⁻ h, (srRun f c N h).1 + (srRun f c N h).2.1 ∂(Measure.pi fun _ : Fin N => μ) :=
        lintegral_mono fun h => (srRun_count f c N h).trans le_add_self
    _ ≤ N := lintegral_srRun_le μ c hf h1 N

/-- CUSUM against Shiryaev–Roberts, step by step: CUSUM has had fewer alarms, or as many with its
statistic no larger. -/
theorem cusum_sr_invariant : ∀ (n : ℕ) (h : Fin n → ℝ),
    (cusumRun f c n h).2 < (srRun f c n h).2.2 ∨
      ((cusumRun f c n h).2 = (srRun f c n h).2.2 ∧ (cusumRun f c n h).1 ≤ (srRun f c n h).1)
  | 0, _ => Or.inr ⟨rfl, le_rfl⟩
  | n + 1, h => by
    set b := f n (Fin.init h) (h (Fin.last n))
    have ih := cusum_sr_invariant n (Fin.init h)
    simp only [cusumRun, srRun]
    set s := srRun f c n (Fin.init h)
    set g := cusumRun f c n (Fin.init h)
    unfold cusumStep srStep
    rcases ih with hlt | ⟨heq, hle⟩
    · split_ifs with h₁ h₂ h₂
      · exact Or.inl (by dsimp only; omega)
      · rcases Nat.lt_or_ge (g.2 + 1) s.2.2 with h₃ | h₃
        · exact Or.inl (by dsimp only; omega)
        · exact Or.inr ⟨by dsimp only; omega, zero_le⟩
      · exact Or.inl (by dsimp only; omega)
      · exact Or.inl (by dsimp only; omega)
    · have hmax : b * max g.1 1 ≤ b * (s.1 + 1) := by
        gcongr
        exact max_le (hle.trans le_self_add) le_add_self
      split_ifs with h₁ h₂ h₂
      · exact Or.inr ⟨by simp only; omega, le_rfl⟩
      · exact absurd (h₁.trans hmax) h₂
      · exact Or.inl (by simp only; omega)
      · exact Or.inr ⟨heq, hmax⟩

/-- CUSUM alarms no more often than Shiryaev–Roberts, on every path. -/
theorem cusum_count_le (n : ℕ) (h : Fin n → ℝ) : (cusumRun f c n h).2 ≤ (srRun f c n h).2.2 := by
  rcases cusum_sr_invariant f c n h with hlt | ⟨heq, -⟩
  · exact hlt.le
  · exact heq.le

/-- One Shiryaev–Roberts step frozen at the first alarm, from `(Y, alarmed, T)`: while running,
`Y = b (Y + 1)`, the first alarm at `c ≤ Y`, and `T` counts the step; once alarmed, nothing moves.
So `T` is `τ ∧ n`, the steps before the first alarm. -/
noncomputable def srFirstStep (b : ℝ≥0∞) (s : ℝ≥0∞ × Bool × ℕ) : ℝ≥0∞ × Bool × ℕ :=
  if s.2.1 then s
  else if c ≤ b * (s.1 + 1) then (b * (s.1 + 1), true, s.2.2 + 1) else (b * (s.1 + 1), false,
      s.2.2 + 1)

/-- Shiryaev–Roberts up to its first alarm. -/
noncomputable def srFirst : (n : ℕ) → (Fin n → ℝ) → ℝ≥0∞ × Bool × ℕ
  | 0, _ => (0, false, 0)
  | n + 1, h => srFirstStep c (f n (Fin.init h) (h (Fin.last n))) (srFirst n (Fin.init h))

/-- One CUSUM step frozen at the first alarm. -/
noncomputable def cusumFirstStep (b : ℝ≥0∞) (s : ℝ≥0∞ × Bool × ℕ) : ℝ≥0∞ × Bool × ℕ :=
  if s.2.1 then s
  else if c ≤ b * max s.1 1 then (b * max s.1 1, true, s.2.2 + 1) else (b * max s.1 1, false,
      s.2.2 + 1)

/-- CUSUM up to its first alarm. -/
noncomputable def cusumFirst : (n : ℕ) → (Fin n → ℝ) → ℝ≥0∞ × Bool × ℕ
  | 0, _ => (0, false, 0)
  | n + 1, h => cusumFirstStep c (f n (Fin.init h) (h (Fin.last n))) (cusumFirst n (Fin.init h))

theorem srFirst_alarm : ∀ (n : ℕ) (h : Fin n → ℝ), (srFirst f c n h).2.1 = true →
    c ≤ (srFirst f c n h).1
  | 0, _, hh => by simp [srFirst] at hh
  | n + 1, h, hh => by
    have ih := srFirst_alarm n (Fin.init h)
    simp only [srFirst] at hh ⊢
    unfold srFirstStep at hh ⊢
    split_ifs at hh ⊢ with h₁ h₂ <;> first | exact ih h₁ | exact h₂

theorem measurable_srFirstStep {α : Type*} [MeasurableSpace α] {b : α → ℝ≥0∞}
    {s : α → ℝ≥0∞ × Bool × ℕ} (hb : Measurable b) (hs : Measurable s) :
    Measurable fun a => srFirstStep c (b a) (s a) := by
  unfold srFirstStep
  have hR : Measurable fun a => b a * ((s a).1 + 1) := hb.mul (hs.fst.add measurable_const)
  have hT : Measurable fun a => (s a).2.2 + 1 :=
    (measurable_from_top (f := fun k : ℕ => k + 1)).comp hs.snd.snd
  refine Measurable.ite ?_ hs (Measurable.ite (measurableSet_le measurable_const hR)
    (hR.prodMk (measurable_const.prodMk hT)) (hR.prodMk (measurable_const.prodMk hT)))
  exact hs.snd.fst (measurableSet_singleton true)

variable {f} in
theorem measurable_srFirst (hf : ∀ n, Measurable (Function.uncurry (f n))) :
    ∀ n, Measurable (srFirst f c n)
  | 0 => measurable_const
  | n + 1 => by
    have hb := (hf n).comp ((measurable_init n).prodMk (measurable_pi_apply (Fin.last n)))
    exact measurable_srFirstStep c hb ((measurable_srFirst hf n).comp (measurable_init n))

/-- One step frozen at the first alarm adds, in expectation, at most the step it counts. -/
theorem srFirst_step_le [IsProbabilityMeasure μ] (hfx : Measurable fun x => f n h x)
    (h1 : ∫⁻ x, f n h x ∂μ ≤ 1) (s : ℝ≥0∞ × Bool × ℕ) :
    ∫⁻ x, (srFirstStep c (f n h x) s).1 ∂μ ≤ s.1 + (if s.2.1 then 0 else 1) := by
  unfold srFirstStep
  by_cases hs : s.2.1 = true
  · simp [hs]
  · have hs' : s.2.1 = false := by simpa using hs
    simp only [hs', Bool.false_eq_true, if_false]
    have : ∀ x, (if c ≤ f n h x * (s.1 + 1) then (f n h x * (s.1 + 1), true, s.2.2 + 1)
        else (f n h x * (s.1 + 1), false, s.2.2 + 1)).1 = f n h x * (s.1 + 1) := fun x => by
      split_ifs <;> rfl
    simp only [this]
    rw [lintegral_mul_const _ hfx]
    calc (∫⁻ x, f n h x ∂μ) * (s.1 + 1) ≤ 1 * (s.1 + 1) := by gcongr
      _ = s.1 + 1 := one_mul _

theorem srFirst_count (b : ℝ≥0∞) (s : ℝ≥0∞ × Bool × ℕ) :
    (((srFirstStep c b s).2.2 : ℕ) : ℝ≥0∞) = s.2.2 + (if s.2.1 then 0 else 1) := by
  unfold srFirstStep
  by_cases hs : s.2.1 = true
  · simp [hs]
  · have hs' : s.2.1 = false := by simpa using hs
    simp only [hs', Bool.false_eq_true, if_false]
    split_ifs <;> push_cast <;> rfl

variable {f} in
/-- **Frozen at its first alarm, Shiryaev–Roberts's statistic is worth at most the steps taken:**
`E[R_(τ∧n)] ≤ E[τ ∧ n]`, `Rₙ − n` being a supermartingale. -/
theorem lintegral_srFirst_le [IsProbabilityMeasure μ] (hf : ∀ n, Measurable (Function.uncurry (f
    n)))
    (h1 : ∀ n h, ∫⁻ x, f n h x ∂μ ≤ 1) :
    ∀ n, ∫⁻ h, (srFirst f c n h).1 ∂(Measure.pi fun _ : Fin n => μ) ≤
      ∫⁻ h, ((srFirst f c n h).2.2 : ℝ≥0∞) ∂(Measure.pi fun _ : Fin n => μ)
  | 0 => by simp [srFirst]
  | n + 1 => by
    have hm := measurable_srFirst c hf
    have hT : ∀ k, Measurable fun h => ((srFirst f c k h).2.2 : ℝ≥0∞) := fun k =>
      (measurable_from_top (f := fun j : ℕ => (j : ℝ≥0∞))).comp (hm k).snd.snd
    set D : (Fin n → ℝ) → ℝ≥0∞ := fun h => if (srFirst f c n h).2.1 then 0 else 1
    have hD : Measurable D := by
      refine Measurable.ite ?_ measurable_const measurable_const
      exact (hm n).snd.fst (measurableSet_singleton true)
    rw [lintegral_pi_succ μ n (g := fun h => (srFirst f c (n + 1) h).1) (hm (n + 1)).fst,
      lintegral_pi_succ μ n (g := fun h => ((srFirst f c (n + 1) h).2.2 : ℝ≥0∞)) (hT (n + 1))]
    have hY : ∀ h : Fin n → ℝ, ∫⁻ x, (srFirst f c (n + 1) (Fin.snoc h x)).1 ∂μ ≤
        (srFirst f c n h).1 + D h := fun h => by
      simp only [srFirst, Fin.init_snoc, Fin.snoc_last]
      exact srFirst_step_le μ f c ((hf n).comp (measurable_const.prodMk measurable_id)) (h1 n h) _
    have hTe : ∀ h : Fin n → ℝ, ∫⁻ x, ((srFirst f c (n + 1) (Fin.snoc h x)).2.2 : ℝ≥0∞) ∂μ =
        (srFirst f c n h).2.2 + D h := fun h => by
      simp only [srFirst, Fin.init_snoc, Fin.snoc_last, srFirst_count]
      rw [lintegral_const, measure_univ, mul_one]
    calc ∫⁻ h, ∫⁻ x, (srFirst f c (n + 1) (Fin.snoc h x)).1 ∂μ ∂(Measure.pi fun _ : Fin n => μ)
        ≤ ∫⁻ h, (srFirst f c n h).1 + D h ∂(Measure.pi fun _ : Fin n => μ) := lintegral_mono hY
      _ = ∫⁻ h, (srFirst f c n h).1 ∂(Measure.pi fun _ : Fin n => μ) +
            ∫⁻ h, D h ∂(Measure.pi fun _ : Fin n => μ) := lintegral_add_right _ hD
      _ ≤ ∫⁻ h, ((srFirst f c n h).2.2 : ℝ≥0∞) ∂(Measure.pi fun _ : Fin n => μ) +
            ∫⁻ h, D h ∂(Measure.pi fun _ : Fin n => μ) := by
          exact add_le_add (lintegral_srFirst_le hf h1 n) le_rfl
      _ = ∫⁻ h, ((srFirst f c n h).2.2 : ℝ≥0∞) + D h ∂(Measure.pi fun _ : Fin n => μ) :=
          (lintegral_add_right _ hD).symm
      _ = ∫⁻ h, ∫⁻ x, ((srFirst f c (n + 1) (Fin.snoc h x)).2.2 : ℝ≥0∞) ∂μ
            ∂(Measure.pi fun _ : Fin n => μ) := by simp only [hTe]

variable {f} in
/-- **Shiryaev–Roberts's mean run length, at a horizon** (Vovk 2021, Proposition 4.4): under the
null,
`c P(τ ≤ N) ≤ E[τ ∧ N]` for its first alarm `τ`. -/
theorem sr_firstAlarm_le [IsProbabilityMeasure μ] (hf : ∀ n, Measurable (Function.uncurry (f n)))
    (h1 : ∀ n h, ∫⁻ x, f n h x ∂μ ≤ 1) (N : ℕ) :
    c * (Measure.pi fun _ : Fin N => μ) {h | (srFirst f c N h).2.1 = true} ≤
      ∫⁻ h, ((srFirst f c N h).2.2 : ℝ≥0∞) ∂(Measure.pi fun _ : Fin N => μ) := by
  have hm := measurable_srFirst c hf N
  have hS : MeasurableSet {h : Fin N → ℝ | (srFirst f c N h).2.1 = true} :=
    hm.snd.fst (measurableSet_singleton true)
  calc c * (Measure.pi fun _ : Fin N => μ) {h | (srFirst f c N h).2.1 = true}
      = ∫⁻ h, Set.indicator {h | (srFirst f c N h).2.1 = true} (fun _ => c) h
          ∂(Measure.pi fun _ : Fin N => μ) := by
        rw [lintegral_indicator_const hS]
    _ ≤ ∫⁻ h, (srFirst f c N h).1 ∂(Measure.pi fun _ : Fin N => μ) := by
        refine lintegral_mono fun h => ?_
        by_cases hh : (srFirst f c N h).2.1 = true
        · rw [Set.indicator_of_mem (show h ∈ {h : Fin N → ℝ | (srFirst f c N h).2.1 = true} from
            hh)]
          exact srFirst_alarm f c N h hh
        · rw [Set.indicator_of_notMem
            (show h ∉ {h : Fin N → ℝ | (srFirst f c N h).2.1 = true} from hh)]
          exact zero_le
    _ ≤ _ := lintegral_srFirst_le μ c hf h1 N

/-- Up to their first alarms, CUSUM's statistic stays under Shiryaev–Roberts's, so CUSUM's first
alarm
comes no sooner. -/
theorem cusumFirst_sr_invariant : ∀ (n : ℕ) (h : Fin n → ℝ),
    (srFirst f c n h).2.2 ≤ (cusumFirst f c n h).2.2 ∧ ((srFirst f c n h).2.1 = true ∨
      ((cusumFirst f c n h).2.1 = false ∧ (cusumFirst f c n h).1 ≤ (srFirst f c n h).1))
  | 0, _ => ⟨le_rfl, Or.inr ⟨rfl, le_rfl⟩⟩
  | n + 1, h => by
    set b := f n (Fin.init h) (h (Fin.last n))
    obtain ⟨hT, hI⟩ := cusumFirst_sr_invariant n (Fin.init h)
    simp only [srFirst, cusumFirst]
    set s := srFirst f c n (Fin.init h)
    set g := cusumFirst f c n (Fin.init h)
    unfold srFirstStep cusumFirstStep
    rcases hI with hA | ⟨hC, hle⟩
    · simp only [hA, if_true]
      refine ⟨?_, Or.inl trivial⟩
      split_ifs <;> (try dsimp only) <;> omega
    · have hmax : b * max g.1 1 ≤ b * (s.1 + 1) := by
        gcongr
        exact max_le (hle.trans le_self_add) le_add_self
      by_cases hsA : s.2.1 = true
      · simp only [hsA, if_true, hC, Bool.false_eq_true, if_false]
        refine ⟨?_, Or.inl trivial⟩
        split_ifs <;> (try dsimp only) <;> omega
      · have hsA' : s.2.1 = false := by simpa using hsA
        simp only [hsA', hC, Bool.false_eq_true, if_false]
        split_ifs with h₁ h₂ h₂
        · exact ⟨by (try dsimp only); omega, Or.inl rfl⟩
        · exact ⟨by (try dsimp only); omega, Or.inl rfl⟩
        · exact absurd (h₂.trans hmax) h₁
        · exact ⟨by (try dsimp only); omega, Or.inr ⟨rfl, hmax⟩⟩

theorem measurable_cusumFirstStep {α : Type*} [MeasurableSpace α] {b : α → ℝ≥0∞}
    {s : α → ℝ≥0∞ × Bool × ℕ} (hb : Measurable b) (hs : Measurable s) :
    Measurable fun a => cusumFirstStep c (b a) (s a) := by
  unfold cusumFirstStep
  have hR : Measurable fun a => b a * max (s a).1 1 := hb.mul (hs.fst.max measurable_const)
  have hT : Measurable fun a => (s a).2.2 + 1 :=
    (measurable_from_top (f := fun k : ℕ => k + 1)).comp hs.snd.snd
  refine Measurable.ite ?_ hs (Measurable.ite (measurableSet_le measurable_const hR)
    (hR.prodMk (measurable_const.prodMk hT)) (hR.prodMk (measurable_const.prodMk hT)))
  exact hs.snd.fst (measurableSet_singleton true)

variable {f} in
theorem measurable_cusumFirst (hf : ∀ n, Measurable (Function.uncurry (f n))) :
    ∀ n, Measurable (cusumFirst f c n)
  | 0 => measurable_const
  | n + 1 => by
    have hb := (hf n).comp ((measurable_init n).prodMk (measurable_pi_apply (Fin.last n)))
    exact measurable_cusumFirstStep c hb ((measurable_cusumFirst hf n).comp (measurable_init n))

variable {f} in
/-- **CUSUM's mean run length, at a horizon** (Vovk 2021, Proposition 4.2): `c P(τ ≤ N) ≤ E[τ ∧ N]`
for its first alarm, which comes no sooner than Shiryaev–Roberts's. -/
theorem cusum_firstAlarm_le [IsProbabilityMeasure μ] (hf : ∀ n, Measurable (Function.uncurry (f n)))
    (h1 : ∀ n h, ∫⁻ x, f n h x ∂μ ≤ 1) (N : ℕ) :
    c * (Measure.pi fun _ : Fin N => μ) {h | (cusumFirst f c N h).2.1 = true} ≤
      ∫⁻ h, ((cusumFirst f c N h).2.2 : ℝ≥0∞) ∂(Measure.pi fun _ : Fin N => μ) := by
  have hsub : {h : Fin N → ℝ | (cusumFirst f c N h).2.1 = true} ⊆
      {h | (srFirst f c N h).2.1 = true} := fun h hh => by
    rcases (cusumFirst_sr_invariant f c N h).2 with hA | ⟨hC, -⟩
    · exact hA
    · simp only [Set.mem_setOf_eq] at hh
      rw [hh] at hC
      exact absurd hC (by decide)
  calc c * (Measure.pi fun _ : Fin N => μ) {h | (cusumFirst f c N h).2.1 = true}
      ≤ c * (Measure.pi fun _ : Fin N => μ) {h | (srFirst f c N h).2.1 = true} := by
        gcongr
    _ ≤ ∫⁻ h, ((srFirst f c N h).2.2 : ℝ≥0∞) ∂(Measure.pi fun _ : Fin N => μ) :=
        sr_firstAlarm_le μ c hf h1 N
    _ ≤ _ := lintegral_mono fun h => Nat.mono_cast (cusumFirst_sr_invariant f c N h).1

variable {f} in
/-- **CUSUM's false alarms** (Vovk 2021, Corollary 4.3, in expectation at every horizon): under
the null, the
expected
number of alarms among the first `N` arrivals is at most `N/c`. -/
theorem cusum_alarms_le [IsProbabilityMeasure μ] (hf : ∀ n, Measurable (Function.uncurry (f n)))
    (h1 : ∀ n h, ∫⁻ x, f n h x ∂μ ≤ 1) (N : ℕ) :
    c * ∫⁻ h, ((cusumRun f c N h).2 : ℝ≥0∞) ∂(Measure.pi fun _ : Fin N => μ) ≤ N :=
  calc c * ∫⁻ h, ((cusumRun f c N h).2 : ℝ≥0∞) ∂(Measure.pi fun _ : Fin N => μ)
      ≤ c * ∫⁻ h, ((srRun f c N h).2.2 : ℝ≥0∞) ∂(Measure.pi fun _ : Fin N => μ) := by
        gcongr with h
        exact_mod_cast cusum_count_le f c N h
    _ ≤ N := sr_alarms_le μ c hf h1 N

/-- Once alarmed, the frozen run stays alarmed. -/
theorem srFirst_alarm_succ (n : ℕ) (h : Fin (n + 1) → ℝ) (hh : (srFirst f c n (Fin.init h)).2.1
    = true) :
    (srFirst f c (n + 1) h).2.1 = true := by
  simp only [srFirst]
  unfold srFirstStep
  simp [hh]

/-- Along a path, the frozen run's step count is the number of horizons before its first alarm. -/
theorem srFirst_count_path (x : ℕ → ℝ) : ∀ N, ((srFirst f c N (fun i : Fin N => x i)).2.2 : ℝ≥0∞) =
    ∑ n ∈ Finset.range N, if (srFirst f c n (fun i : Fin n => x i)).2.1 then 0 else 1
  | 0 => by simp [srFirst]
  | N + 1 => by
    have hinit : Fin.init (fun i : Fin (N + 1) => x i) = fun i : Fin N => x i := rfl
    rw [Finset.sum_range_succ, ← srFirst_count_path x N]
    simp only [srFirst, hinit]
    exact srFirst_count c _ _

/-- The measure fact behind a mean run length: alarm events `E n`, increasing, with
`c P(E N) ≤ ∑_{n < N} P(not yet)` at every horizon, give `c ≤ ∑ₙ P(not yet at n)`, the expected
time to the alarm. -/
theorem le_tsum_of_alarm_bound {Ω : Type*} [MeasurableSpace Ω] (P : Measure Ω)
    [IsProbabilityMeasure P]
    (E : ℕ → Set Ω) (hE : Monotone E) (hm : ∀ n, MeasurableSet (E n))
    (hb : ∀ N, c * P (E N) ≤ ∑ n ∈ Finset.range N, P (E n)ᶜ) : c ≤ ∑' n, P (E n)ᶜ := by
  by_cases h0 : P (⋃ n, E n)ᶜ = 0
  · have hU : P (⋃ n, E n) = 1 := by
      rw [prob_compl_eq_zero_iff (MeasurableSet.iUnion hm)] at h0
      exact h0
    have ht := ENNReal.Tendsto.const_mul (tendsto_measure_iUnion_atTop hE) (a := c)
      (Or.inl (by rw [hU]; exact one_ne_zero))
    rw [hU, mul_one] at ht
    exact le_of_tendsto' ht fun N => (hb N).trans (ENNReal.sum_le_tsum _)
  · have hle : ∀ n, P (⋃ n, E n)ᶜ ≤ P (E n)ᶜ := fun n =>
      measure_mono (Set.compl_subset_compl.2 (Set.subset_iUnion E n))
    calc c ≤ ∞ := le_top
      _ = ∑' _ : ℕ, P (⋃ n, E n)ᶜ := (ENNReal.tsum_const_eq_top_of_ne_zero h0).symm
      _ ≤ ∑' n, P (E n)ᶜ := ENNReal.tsum_le_tsum hle

variable {f} in
open Filter Topology in
/-- **Shiryaev–Roberts's mean time to a false alarm is at least `c`** (Vovk 2021, Proposition 4.4).
Let the p-values `p i` on `Ω` be independent with law `μ` at every horizon, and each bet be
measurable, integrate to at most 1 and read only the earlier p-values. Then the expected time to the
first alarm, `∑ₙ 1{no alarm by n}` (infinite if it never alarms), is at least `c`. -/
theorem sr_meanRunLength [IsProbabilityMeasure μ] (hf : ∀ n, Measurable (Function.uncurry (f n)))
    (h1 : ∀ n h, ∫⁻ x, f n h x ∂μ ≤ 1) {Ω : Type*} [MeasurableSpace Ω] (P : Measure Ω)
    [IsProbabilityMeasure P] (p : ℕ → Ω → ℝ) (hp : ∀ i, Measurable (p i))
    (hlaw : ∀ n, P.map (fun ω (i : Fin n) => p i ω) = Measure.pi fun _ => μ) :
    c ≤ ∫⁻ ω, ∑' n, (if (srFirst f c n (fun i : Fin n => p i ω)).2.1 then 0 else 1) ∂P := by
  have hφ : ∀ n, Measurable fun ω (i : Fin n) => p i ω := fun n =>
    measurable_pi_lambda _ fun i => hp i
  set E : ℕ → Set Ω := fun n => {ω | (srFirst f c n (fun i : Fin n => p i ω)).2.1 = true}
  have hmE : ∀ n, MeasurableSet (E n) := fun n =>
    ((measurable_srFirst c hf n).comp (hφ n)).snd.fst (measurableSet_singleton true)
  have hind : ∀ n ω, (if (srFirst f c n (fun i : Fin n => p i ω)).2.1 then (0 : ℝ≥0∞) else 1) =
      Set.indicator (E n)ᶜ (fun _ => 1) ω := fun n ω => by
    by_cases hω : (srFirst f c n (fun i : Fin n => p i ω)).2.1 = true
    · rw [if_pos hω, Set.indicator_of_notMem (by simpa [E] using hω)]
    · rw [if_neg hω, Set.indicator_of_mem (by simpa [E] using hω)]
  have hint : ∀ n, ∫⁻ ω, Set.indicator (E n)ᶜ (fun _ => (1 : ℝ≥0∞)) ω ∂P = P (E n)ᶜ := fun n => by
    rw [lintegral_indicator_const (hmE n).compl, one_mul]
  simp only [hind]
  rw [lintegral_tsum fun n => (measurable_const.indicator (hmE n).compl).aemeasurable]
  simp only [hint]
  refine le_tsum_of_alarm_bound c P E ?_ hmE fun N => ?_
  · refine monotone_nat_of_le_succ fun n ω hω => ?_
    exact srFirst_alarm_succ f c n (fun i : Fin (n + 1) => p i ω) hω
  · have hmap : ∀ {g : (Fin N → ℝ) → ℝ≥0∞}, Measurable g →
        ∫⁻ ω, g (fun i : Fin N => p i ω) ∂P = ∫⁻ h, g h ∂(Measure.pi fun _ : Fin N => μ) :=
      fun hg => by rw [← hlaw N, lintegral_map hg (hφ N)]
    have hPE : P (E N) = (Measure.pi fun _ : Fin N => μ) {h | (srFirst f c N h).2.1 = true} := by
      rw [← hlaw N, Measure.map_apply (hφ N)]
      · rfl
      · exact (measurable_srFirst c hf N).snd.fst (measurableSet_singleton true)
    have hT : ∫⁻ h, ((srFirst f c N h).2.2 : ℝ≥0∞) ∂(Measure.pi fun _ : Fin N => μ) =
        ∑ n ∈ Finset.range N, P (E n)ᶜ := by
      rw [← hmap (g := fun h => ((srFirst f c N h).2.2 : ℝ≥0∞))
        ((measurable_from_top (f := fun j : ℕ => (j : ℝ≥0∞))).comp (measurable_srFirst c hf
            N).snd.snd)]
      have hpt : ∀ ω, ((srFirst f c N (fun i : Fin N => p i ω)).2.2 : ℝ≥0∞) =
          ∑ n ∈ Finset.range N, if (srFirst f c n (fun i : Fin n => p i ω)).2.1 then 0 else 1 :=
        fun ω => srFirst_count_path f c (fun i => p i ω) N
      simp only [hpt]
      simp only [hind]
      rw [lintegral_finsetSum _ fun n _ => measurable_const.indicator (hmE n).compl]
      simp only [hint]
    rw [hPE, ← hT]
    exact sr_firstAlarm_le μ c hf h1 N

variable {f} in
/-- **CUSUM's mean time to a false alarm is at least `c`** (Vovk 2021, Proposition 4.2): its first
alarm comes no sooner than Shiryaev–Roberts's. -/
theorem cusum_meanRunLength [IsProbabilityMeasure μ] (hf : ∀ n, Measurable (Function.uncurry (f n)))
    (h1 : ∀ n h, ∫⁻ x, f n h x ∂μ ≤ 1) {Ω : Type*} [MeasurableSpace Ω] (P : Measure Ω)
    [IsProbabilityMeasure P] (p : ℕ → Ω → ℝ) (hp : ∀ i, Measurable (p i))
    (hlaw : ∀ n, P.map (fun ω (i : Fin n) => p i ω) = Measure.pi fun _ => μ) :
    c ≤ ∫⁻ ω, ∑' n, (if (cusumFirst f c n (fun i : Fin n => p i ω)).2.1 then 0 else 1) ∂P := by
  refine (sr_meanRunLength μ c hf h1 P p hp hlaw).trans (lintegral_mono fun ω =>
    ENNReal.tsum_le_tsum fun n => ?_)
  by_cases hs : (srFirst f c n (fun i : Fin n => p i ω)).2.1 = true
  · simp [hs]
  · rcases (cusumFirst_sr_invariant f c n (fun i : Fin n => p i ω)).2 with hA | ⟨hC, -⟩
    · exact absurd hA hs
    · simp [hs, hC]

variable {f} in
/-- **The false alarms on a stream** (Vovk 2021, Corollaries 4.3 and 4.6, in expectation at
every horizon):
with the
p-values independent with law `μ`, the expected numbers of Shiryaev–Roberts's and CUSUM's alarms
among the first `N` arrivals are at most `N/c`. -/
theorem alarms_stream_le [IsProbabilityMeasure μ] (hf : ∀ n, Measurable (Function.uncurry (f n)))
    (h1 : ∀ n h, ∫⁻ x, f n h x ∂μ ≤ 1) {Ω : Type*} [MeasurableSpace Ω] (P : Measure Ω)
    (p : ℕ → Ω → ℝ) (hp : ∀ i, Measurable (p i))
    (hlaw : ∀ n, P.map (fun ω (i : Fin n) => p i ω) = Measure.pi fun _ => μ) (N : ℕ) :
    c * ∫⁻ ω, ((cusumRun f c N (fun i : Fin N => p i ω)).2 : ℝ≥0∞) ∂P ≤
        c * ∫⁻ ω, ((srRun f c N (fun i : Fin N => p i ω)).2.2 : ℝ≥0∞) ∂P ∧
      c * ∫⁻ ω, ((srRun f c N (fun i : Fin N => p i ω)).2.2 : ℝ≥0∞) ∂P ≤ N := by
  have hφ : Measurable fun ω (i : Fin N) => p i ω := measurable_pi_lambda _ fun i => hp i
  have hA : Measurable fun h => ((srRun f c N h).2.2 : ℝ≥0∞) :=
    (measurable_from_top (f := fun k : ℕ => (k : ℝ≥0∞))).comp (measurable_srRun c hf N).snd.snd
  refine ⟨?_, ?_⟩
  · gcongr with ω
    exact Nat.mono_cast (cusum_count_le f c N _)
  · have : ∫⁻ ω, ((srRun f c N (fun i : Fin N => p i ω)).2.2 : ℝ≥0∞) ∂P =
        ∫⁻ h, ((srRun f c N h).2.2 : ℝ≥0∞) ∂(Measure.pi fun _ : Fin N => μ) := by
      rw [← hlaw N, lintegral_map hA hφ]
    rw [this]
    exact sr_alarms_le μ c hf h1 N

end ChangeAlarms

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

/-! ### Slice 10: the full step's error, completed -/

section FullStepError

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

omit [Fintype ι] [DecidableEq ι] in
/-- The share of units without a property is one less the share with it. -/
theorem card_filter_not_div {C : Finset ι} (hC : C.Nonempty) (q : ι → Prop) [DecidablePred q] :
    ((C.filter fun u => ¬ q u).card : ℝ) / C.card = 1 - ((C.filter q).card : ℝ) / C.card := by
  have hc : (0 : ℝ) < C.card := by exact_mod_cast card_pos.2 hC
  have h : ((C.filter q).card : ℝ) + ((C.filter fun u => ¬ q u).card : ℝ) = C.card := by
    exact_mod_cast C.card_filter_add_card_filter_not q
  rw [eq_sub_iff_add_eq, ← add_div, add_comm, h, div_self hc.ne']

omit [Fintype ι] [DecidableEq ι] in
/-- Its weighted form, over draws with weights summing to one. -/
theorem weighted_share_not {W : Type*} [Fintype W] {C : Finset ι} (hC : C.Nonempty)
    (p : W → ℝ) (hp1 : ∑ w, p w = 1) (bad : W → ι → Prop) [∀ w, DecidablePred (bad w)] :
    ∑ w, p w * (((C.filter fun u => ¬ bad w u).card : ℝ) / C.card) =
      1 - ∑ w, p w * (((C.filter (bad w)).card : ℝ) / C.card) := by
  simp only [card_filter_not_div hC, mul_sub, mul_one, sum_sub_distrib, hp1]

/-- `ν` as a share: a weighted count of orderings at most `ν` times their number is a weighted
share at most `ν`. -/
theorem weighted_share_le {W : Type*} [Fintype W] (C : Finset ι) (p : W → ℝ)
    (bad : W → Equiv.Perm ι → Prop) [∀ w, DecidablePred (bad w)] {ν : ℝ}
    (hfail : ∑ w, p w * (((permsOn C).filter (bad w)).card : ℝ) ≤ ν * (permsOn C).card) :
    ∑ w, p w * ((((permsOn C).filter (bad w)).card : ℝ) / (permsOn C).card) ≤ ν := by
  have hg : (0 : ℝ) < (permsOn C).card := by exact_mod_cast card_permsOn_pos C
  have heq : ∑ w, p w * ((((permsOn C).filter (bad w)).card : ℝ) / (permsOn C).card) =
      (∑ w, p w * (((permsOn C).filter (bad w)).card : ℝ)) / (permsOn C).card := by
    rw [sum_div]
    exact sum_congr rfl fun w _ => by ring
  rw [heq, div_le_iff₀ hg]
  exact hfail

/-- Fewer orderings fail less. -/
theorem weighted_fail_mono {W : Type*} [Fintype W] (C : Finset ι) (p : W → ℝ)
    (hp0 : ∀ w, 0 ≤ p w) (a b : W → Equiv.Perm ι → Prop) [∀ w, DecidablePred (a w)]
    [∀ w, DecidablePred (b w)] (hab : ∀ w σ, a w σ → b w σ) {ν : ℝ}
    (h : ∑ w, p w * (((permsOn C).filter (b w)).card : ℝ) ≤ ν * (permsOn C).card) :
    ∑ w, p w * (((permsOn C).filter (a w)).card : ℝ) ≤ ν * (permsOn C).card := by
  refine le_trans (sum_le_sum fun w _ => mul_le_mul_of_nonneg_left ?_ (hp0 w)) h
  exact_mod_cast card_le_card fun σ hσ => by
    rw [mem_filter] at hσ ⊢
    exact ⟨hσ.1, hab w σ hσ.2⟩

/-- **Coverage from a sufficient condition.** If every ordering that puts the test unit on a unit
with property `good` is accepted, the share of orderings accepted is at least the share of units
with it. -/
theorem coverage_ge_of_good {C : Finset ι} {t : ι} (ht : t ∈ C) (good : ι → Prop)
    [DecidablePred good] (acc : Equiv.Perm ι → Prop) [DecidablePred acc]
    (h : ∀ σ ∈ permsOn C, good (σ t) → acc σ) :
    ((C.filter good).card : ℝ) / C.card ≤
      (((permsOn C).filter acc).card : ℝ) / (permsOn C).card := by
  rw [← ratio_permsOn ht good]
  have hg : (0 : ℝ) < (permsOn C).card := by exact_mod_cast card_permsOn_pos C
  refine div_le_div_of_nonneg_right ?_ hg.le
  exact_mod_cast card_le_card fun σ hσ => by
    rw [mem_filter] at hσ ⊢
    exact ⟨hσ.1, h σ hσ.1 hσ.2⟩

/-- **The same, with a chance of failure.** `W` carries the draws' randomness, with weights `p`.
The draw `w` is trusted on the orderings `hold w σ`, and on the others its estimates may be
anything. If every trusted ordering that puts the test unit on a good unit is accepted, coverage is
at least the good units' share, less `ν`, the weighted share of untrusted orderings among those. -/
theorem coverage_ge_of_good_fail {W : Type*} [Fintype W] {C : Finset ι} {t : ι} (ht : t ∈ C)
    (p : W → ℝ) (hp0 : ∀ w, 0 ≤ p w) (good : W → ι → Prop) [∀ w, DecidablePred (good w)]
    (hold acc : W → Equiv.Perm ι → Prop) [∀ w, DecidablePred (hold w)]
    [∀ w, DecidablePred (acc w)] (h : ∀ w, ∀ σ ∈ permsOn C, good w (σ t) → hold w σ → acc w σ)
    {ν : ℝ} (hfail : ∑ w, p w * (((permsOn C).filter fun σ : Equiv.Perm ι =>
        good w (σ t) ∧ ¬ hold w σ).card : ℝ) ≤ ν * (permsOn C).card) :
    ∑ w, p w * (((C.filter (good w)).card : ℝ) / C.card) - ν ≤
      ∑ w, p w * ((((permsOn C).filter (acc w)).card : ℝ) / (permsOn C).card) := by
  have hg : (0 : ℝ) < (permsOn C).card := by exact_mod_cast card_permsOn_pos C
  have hw : ∀ w, ((C.filter (good w)).card : ℝ) / C.card ≤
      (((permsOn C).filter (acc w)).card : ℝ) / (permsOn C).card +
        (((permsOn C).filter fun σ : Equiv.Perm ι => good w (σ t) ∧ ¬ hold w σ).card : ℝ) /
          (permsOn C).card := by
    intro w
    rw [← ratio_permsOn ht (good w), ← add_div]
    refine div_le_div_of_nonneg_right ?_ hg.le
    have hsub : ((permsOn C).filter fun σ : Equiv.Perm ι => good w (σ t)) ⊆
        ((permsOn C).filter (acc w)) ∪
          (permsOn C).filter fun σ : Equiv.Perm ι => good w (σ t) ∧ ¬ hold w σ := by
      intro σ hσ
      rw [mem_filter] at hσ
      rw [mem_union, mem_filter, mem_filter]
      by_cases hh : hold w σ
      · exact Or.inl ⟨hσ.1, h w σ hσ.1 hσ.2 hh⟩
      · exact Or.inr ⟨hσ.1, hσ.2, hh⟩
    exact_mod_cast (card_le_card hsub).trans (card_union_le _ _)
  have hsum : ∑ w, p w * (((C.filter (good w)).card : ℝ) / C.card) ≤
      ∑ w, p w * ((((permsOn C).filter (acc w)).card : ℝ) / (permsOn C).card) +
        ∑ w, p w * ((((permsOn C).filter fun σ : Equiv.Perm ι =>
          good w (σ t) ∧ ¬ hold w σ).card : ℝ) / (permsOn C).card) := by
    rw [← sum_add_distrib]
    refine sum_le_sum fun w _ => ?_
    rw [← mul_add]
    exact mul_le_mul_of_nonneg_left (hw w) (hp0 w)
  have hν : ∑ w, p w * ((((permsOn C).filter fun σ : Equiv.Perm ι =>
      good w (σ t) ∧ ¬ hold w σ).card : ℝ) / (permsOn C).card) ≤ ν :=
    weighted_share_le C p (fun w σ => good w (σ t) ∧ ¬ hold w σ) hfail
  linarith

/-! #### The equi-tailed set -/

/-- **The equi-tailed set at `k` per tail**: CB-1's rule on the score (the lower tail) and on its
negation (the upper). The full step's CDF score `F̃` takes `k = ⌊α·|C|/2⌋`. -/
def AcceptsTwoSided (k : ℕ) (C : Finset ι) (s : ι → ℝ) (u : ι) : Prop :=
  Accepts k C s u ∧ Accepts k C (fun j => -s j) u

noncomputable instance (k : ℕ) (C : Finset ι) (s : ι → ℝ) (u : ι) :
    Decidable (AcceptsTwoSided k C s u) :=
  inferInstanceAs (Decidable (_ ∧ _))

omit [Fintype ι] [DecidableEq ι] in
/-- The exact set refuses at most `2k` units, `k` per tail. -/
theorem card_not_acceptsTwoSided_le (C : Finset ι) (s : ι → ℝ) (k : ℕ) :
    (C.filter fun u => ¬ AcceptsTwoSided k C s u).card ≤ 2 * k := by
  classical
  have hsub : (C.filter fun u => ¬ AcceptsTwoSided k C s u) ⊆
      (C.filter fun u => rankLe C s u ≤ k) ∪
        C.filter fun u => rankLe C (fun j => -s j) u ≤ k := by
    intro u hu
    rw [mem_filter] at hu
    rw [mem_union, mem_filter, mem_filter]
    by_cases h : Accepts k C s u
    · exact Or.inr ⟨hu.1, not_lt.1 fun h' => hu.2 ⟨h, h'⟩⟩
    · exact Or.inl ⟨hu.1, not_lt.1 h⟩
  have h₁ := card_rankLe_le C s k
  have h₂ := card_rankLe_le C (fun j => -s j) k
  have := (card_le_card hsub).trans (card_union_le _ _)
  omega

theorem acceptsTwoSided_comp {k : ℕ} {C : Finset ι} {σ : Equiv.Perm ι} (hσ : σ ∈ permsOn C)
    (s : ι → ℝ) (u : ι) : AcceptsTwoSided k C (s ∘ σ) u ↔ AcceptsTwoSided k C s (σ u) := by
  have h₁ := rankLe_comp hσ s u
  have h₂ := rankLe_comp hσ (fun j => -s j) u
  simp only [AcceptsTwoSided, Accepts]
  rw [h₁, show (fun j => -(s ∘ σ) j) = (fun j => -s j) ∘ σ from rfl, h₂]

/-- The two tails' fragile units, the margin of the equi-tailed set. -/
def FragileTwoSided (k : ℕ) (C : Finset ι) (s : ι → ℝ) (ε : ι → ℝ) (u : ι) : Prop :=
  Fragile k C s ε u ∨ Fragile k C (fun j => -s j) ε u

noncomputable instance (k : ℕ) (C : Finset ι) (s ε : ι → ℝ) (u : ι) :
    Decidable (FragileTwoSided k C s ε u) :=
  inferInstanceAs (Decidable (_ ∨ _))

theorem fragileTwoSided_comp {k : ℕ} {C : Finset ι} {σ : Equiv.Perm ι} (hσ : σ ∈ permsOn C)
    (s ε : ι → ℝ) (u : ι) :
    FragileTwoSided k C (s ∘ σ) (ε ∘ σ) u ↔ FragileTwoSided k C s ε (σ u) := by
  simp only [FragileTwoSided]
  rw [fragile_comp hσ s ε u, show (fun j => -(s ∘ σ) j) = (fun j => -s j) ∘ σ from rfl,
    fragile_comp hσ (fun j => -s j) ε u]

omit [Fintype ι] [DecidableEq ι] in
/-- Bounds on a score's estimates are bounds on its negation's. -/
theorem close_neg {C : Finset ι} {s ŝ ε : ι → ℝ} (h : ∀ j ∈ C, |ŝ j - s j| ≤ ε j) :
    ∀ j ∈ C, |(fun j => -ŝ j) j - (fun j => -s j) j| ≤ ε j := fun j hj => by
  change |-ŝ j - -s j| ≤ ε j
  rw [show -ŝ j - -s j = -(ŝ j - s j) by ring, abs_neg]
  exact h j hj

omit [Fintype ι] in
/-- **The deterministic part, per tail.** An estimate within its bounds refuses a unit the exact
equi-tailed set accepts only when the unit is fragile in one of the tails. -/
theorem acceptsTwoSided_or_fragile {k : ℕ} {C : Finset ι} {s ŝ ε : ι → ℝ}
    (hclose : ∀ j ∈ C, |ŝ j - s j| ≤ ε j) {u : ι} (hu : u ∈ C) (h : AcceptsTwoSided k C s u) :
    AcceptsTwoSided k C ŝ u ∨ FragileTwoSided k C s ε u := by
  have hclose' := close_neg hclose
  rcases accepts_or_fragile hclose hu h.1 with h₁ | h₁
  · rcases accepts_or_fragile hclose' hu h.2 with h₂ | h₂
    · exact Or.inl ⟨h₁, h₂⟩
    · exact Or.inr (Or.inr h₂)
  · exact Or.inr (Or.inl h₁)

/-- **CB-3 for the equi-tailed set.** With bounds carried by the units, the estimate's equi-tailed
set covers at least `1 − |R ∪ M_lo ∪ M_hi|/|C|`, with `R` the units the exact set refuses (at most
`2k`, `card_not_acceptsTwoSided_le`) and `M_lo`, `M_hi` the two tails' fragile units. The tails'
bands enter as their union, which is never more than their sum. -/
theorem coverage_close_twoSided {Z : Type*} {C : Finset ι} {t : ι} (ht : t ∈ C)
    {S : (ι → Z) → ι → ℝ} (hS : Equivariant C S) (z : ι → Z) (ŝ : Equiv.Perm ι → ι → ℝ)
    (ε : ι → ℝ) (hclose : ∀ σ ∈ permsOn C, ∀ j ∈ C, |ŝ σ j - S (z ∘ σ) j| ≤ ε (σ j)) (k : ℕ) :
    1 - ((C.filter fun u => ¬ AcceptsTwoSided k C (S z) u ∨
        FragileTwoSided k C (S z) ε u).card : ℝ) / C.card ≤
      (((permsOn C).filter fun σ : Equiv.Perm ι => AcceptsTwoSided k C (ŝ σ) t).card : ℝ) /
        (permsOn C).card := by
  classical
  have hgood := coverage_ge_of_good ht
    (fun u => ¬ (¬ AcceptsTwoSided k C (S z) u ∨ FragileTwoSided k C (S z) ε u))
    (fun σ => AcceptsTwoSided k C (ŝ σ) t) fun σ hσ hgood => by
      simp only [not_or, not_not] at hgood
      have hfun : S (z ∘ σ) = S z ∘ σ := funext fun i => hS σ hσ z i
      have hacc : AcceptsTwoSided k C (S (z ∘ σ)) t := by
        rw [hfun, acceptsTwoSided_comp hσ]
        exact hgood.1
      rcases acceptsTwoSided_or_fragile (ε := fun j => ε (σ j)) (hclose σ hσ) ht hacc with h | h
      · exact h
      · exfalso
        apply hgood.2
        rw [hfun] at h
        exact (fragileTwoSided_comp hσ (S z) ε t).1 h
  rw [card_filter_not_div ⟨t, ht⟩] at hgood
  exact hgood

/-- **The inflated equi-tailed set** keeps `1 − 2k/|C|`, with bounds by ordering and position. -/
theorem coverage_inflated_twoSided {Z : Type*} {C : Finset ι} {t : ι} (ht : t ∈ C)
    {S : (ι → Z) → ι → ℝ} (hS : Equivariant C S) (z : ι → Z) (ŝ ε : Equiv.Perm ι → ι → ℝ)
    (hclose : ∀ σ ∈ permsOn C, ∀ j ∈ C, |ŝ σ j - S (z ∘ σ) j| ≤ ε σ j) (k : ℕ) :
    1 - ((2 * k : ℕ) : ℝ) / C.card ≤
      (((permsOn C).filter fun σ : Equiv.Perm ι => AcceptsInflated k C (ŝ σ) (ε σ) t ∧
        AcceptsInflated k C (fun j => -ŝ σ j) (ε σ) t).card : ℝ) / (permsOn C).card := by
  classical
  have hc : (0 : ℝ) < C.card := by exact_mod_cast card_pos.2 ⟨t, ht⟩
  have hgood := coverage_ge_of_good ht (fun u => AcceptsTwoSided k C (S z) u)
    (fun σ => AcceptsInflated k C (ŝ σ) (ε σ) t ∧
      AcceptsInflated k C (fun j => -ŝ σ j) (ε σ) t) fun σ hσ hgood => by
      have hfun : S (z ∘ σ) = S z ∘ σ := funext fun i => hS σ hσ z i
      have hacc : AcceptsTwoSided k C (S (z ∘ σ)) t := by
        rw [hfun, acceptsTwoSided_comp hσ]
        exact hgood
      have hclose' := close_neg (hclose σ hσ)
      exact ⟨acceptsInflated_of_accepts (hclose σ hσ) ht hacc.1,
        acceptsInflated_of_accepts hclose' ht hacc.2⟩
  have hrej := card_not_acceptsTwoSided_le C (S z) k
  have hsplit := C.card_filter_add_card_filter_not (fun u => AcceptsTwoSided k C (S z) u)
  have hcnt : (C.card : ℝ) - ((2 * k : ℕ) : ℝ) ≤ ((C.filter fun u =>
      AcceptsTwoSided k C (S z) u).card : ℝ) := by
    have : C.card ≤ (C.filter fun u => AcceptsTwoSided k C (S z) u).card + 2 * k := by omega
    have : (C.card : ℝ) ≤ ((C.filter fun u => AcceptsTwoSided k C (S z) u).card : ℝ) +
        ((2 * k : ℕ) : ℝ) := by exact_mod_cast this
    linarith
  rw [one_sub_div hc.ne']
  exact (div_le_div_of_nonneg_right hcnt hc.le).trans hgood

omit [Fintype ι] in
/-- **When the tails' bands share a unit.** A unit fragile in both tails has at least
`|C| − 2k + 2` units of `C`, itself among them, within the sum of their bounds of it. At
`k = ⌊α|C|/2⌋` that is most of `C`, so the bands add except where the scores bunch within the
bounds. -/
theorem card_near_of_fragile_both {k : ℕ} {C : Finset ι} {s ε : ι → ℝ} {u : ι} (hu : u ∈ C)
    (hlo : Fragile k C s ε u) (hhi : Fragile k C (fun j => -s j) ε u) :
    C.card + 2 ≤ 2 * k + (C.filter fun j => j = u ∨ |s j - s u| < ε j + ε u).card := by
  classical
  set A := C.filter fun j => j = u ∨ s j + (ε j + ε u) ≤ s u with hA
  set B := C.filter fun j => j = u ∨ -s j + (ε j + ε u) ≤ -s u with hB
  set Nr := C.filter fun j => j = u ∨ |s j - s u| < ε j + ε u with hN
  have hAk : A.card ≤ k := hlo.2
  have hBk : B.card ≤ k := hhi.2
  have huA : u ∈ A := mem_filter.2 ⟨hu, Or.inl rfl⟩
  have huB : u ∈ B := mem_filter.2 ⟨hu, Or.inl rfl⟩
  have hcover : C ⊆ (A.erase u ∪ B.erase u) ∪ Nr := by
    intro j hj
    rw [mem_union, mem_union, mem_erase, mem_erase]
    by_cases hju : j = u
    · exact Or.inr (mem_filter.2 ⟨hj, Or.inl hju⟩)
    · by_cases hnear : |s j - s u| < ε j + ε u
      · exact Or.inr (mem_filter.2 ⟨hj, Or.inr hnear⟩)
      · rw [not_lt, le_abs] at hnear
        rcases hnear with h | h
        · exact Or.inl (Or.inr ⟨hju, mem_filter.2 ⟨hj, Or.inr (by linarith)⟩⟩)
        · exact Or.inl (Or.inl ⟨hju, mem_filter.2 ⟨hj, Or.inr (by linarith)⟩⟩)
  have h1 := card_le_card hcover
  have h2 := card_union_le (A.erase u ∪ B.erase u) Nr
  have h3 := card_union_le (A.erase u) (B.erase u)
  rw [card_erase_of_mem huA, card_erase_of_mem huB] at h3
  have hA1 : 1 ≤ A.card := card_pos.2 ⟨u, huA⟩
  have hB1 : 1 ≤ B.card := card_pos.2 ⟨u, huB⟩
  omega

/-- **The engine's reading, two-sided.** For every ordering, the exact union band is at most the
union band read from that ordering's estimates and the units' bounds. -/
theorem card_fragileTwoSided_le_estimated {Z : Type*} {C : Finset ι} {S : (ι → Z) → ι → ℝ}
    (hS : Equivariant C S) (z : ι → Z) {σ : Equiv.Perm ι} (hσ : σ ∈ permsOn C) (ŝ ε : ι → ℝ)
    (hclose : ∀ j ∈ C, |ŝ j - S (z ∘ σ) j| ≤ ε (σ j)) (k : ℕ) :
    (C.filter fun u => FragileTwoSided k C (S z) ε u).card ≤
      (C.filter fun u => FragileEstimated k C ŝ (ε ∘ σ) u ∨
        FragileEstimated k C (fun j => -ŝ j) (ε ∘ σ) u).card := by
  have heq : (C.filter fun u => FragileTwoSided k C (S z) ε u).card =
      (C.filter fun u => FragileTwoSided k C (S (z ∘ σ)) (ε ∘ σ) u).card := by
    symm
    apply card_equiv σ
    intro u
    have : S (z ∘ σ) = S z ∘ σ := funext fun i => hS σ hσ z i
    simp only [mem_filter, this, fragileTwoSided_comp hσ, mem_iff_of_mem_permsOn hσ]
  rw [heq]
  have hclose' := close_neg (ε := ε ∘ σ) hclose
  refine card_le_card fun u hu => ?_
  rw [mem_filter] at hu ⊢
  refine ⟨hu.1, ?_⟩
  rcases hu.2 with h | h
  · exact Or.inl (fragileEstimated_of_fragile hclose hu.1 h)
  · exact Or.inr (fragileEstimated_of_fragile hclose' hu.1 h)

/-! #### The PIT -/

/-- `P(PIT ≤ k/|C|)` given the scores: the PIT is `(R − 1 + V)/|C|` with `V` uniform and
`R = L + ⌈θ(T + 1)⌉` the randomized rank (`tieRankProb`), so this is `P(R ≤ k)`. -/
noncomputable def pitRankLe (C : Finset ι) (s : ι → ℝ) (u : ι) (k : ℕ) : ℝ :=
  ∑ r ∈ Finset.Icc 1 k, tieRankProb C s u r

omit [Fintype ι] [DecidableEq ι] in
theorem pitRankLe_le (C : Finset ι) (s : ι → ℝ) {u : ι} (hu : u ∈ C) (k : ℕ) :
    pitRankLe C s u k ≤ if rankLt C s u < k then 1 else 0 := by
  classical
  have hlt := rankLt_lt_rankLe C s hu
  have hpos : (0 : ℝ) < ((rankLe C s u - rankLt C s u : ℕ) : ℝ) := by
    exact_mod_cast Nat.sub_pos_of_lt hlt
  unfold pitRankLe tieRankProb
  rw [← sum_filter]
  split_ifs with hk
  · rw [sum_const, nsmul_eq_mul]
    have hsub : ((Finset.Icc 1 k).filter fun r => rankLt C s u < r ∧ r ≤ rankLe C s u) ⊆
        Finset.Ioc (rankLt C s u) (rankLe C s u) := fun r hr => by
      rw [mem_filter] at hr
      exact Finset.mem_Ioc.2 hr.2
    have hcard := card_le_card hsub
    rw [Nat.card_Ioc] at hcard
    rw [mul_one_div, div_le_one hpos]
    exact_mod_cast hcard
  · rw [Finset.sum_eq_zero]
    intro r hr
    rw [mem_filter, Finset.mem_Icc] at hr
    omega

omit [Fintype ι] [DecidableEq ι] in
theorem pitRankLe_nonneg (C : Finset ι) (s : ι → ℝ) (u : ι) (k : ℕ) : 0 ≤ pitRankLe C s u k :=
  sum_nonneg fun r _ => by
    unfold tieRankProb
    split_ifs
    · positivity
    · exact le_rfl

omit [Fintype ι] [DecidableEq ι] in
theorem le_pitRankLe (C : Finset ι) (s : ι → ℝ) {u : ι} (hu : u ∈ C) (k : ℕ) :
    (if rankLe C s u ≤ k then 1 else 0) ≤ pitRankLe C s u k := by
  classical
  have hlt := rankLt_lt_rankLe C s hu
  have hpos : (0 : ℝ) < ((rankLe C s u - rankLt C s u : ℕ) : ℝ) := by
    exact_mod_cast Nat.sub_pos_of_lt hlt
  split_ifs with hk
  · unfold pitRankLe tieRankProb
    rw [← sum_filter]
    have heq : ((Finset.Icc 1 k).filter fun r => rankLt C s u < r ∧ r ≤ rankLe C s u) =
        Finset.Ioc (rankLt C s u) (rankLe C s u) := by
      ext r
      rw [mem_filter, Finset.mem_Icc, Finset.mem_Ioc]
      omega
    rw [heq, sum_const, Nat.card_Ioc, nsmul_eq_mul, mul_one_div, div_self hpos.ne']
  · exact pitRankLe_nonneg C s u k

/-- The units of `C` at least `ε j + ε u` below `u`. -/
noncomputable def shiftedBelow (C : Finset ι) (s ε : ι → ℝ) (u : ι) : ℕ :=
  (C.filter fun j => s j + (ε j + ε u) < s u).card

/-- The units of `C` at most `ε j + ε u` above `u`, `u` among them for bounds `≥ 0`. -/
noncomputable def shiftedAtMost (C : Finset ι) (s ε : ι → ℝ) (u : ι) : ℕ :=
  (C.filter fun j => s j ≤ s u + (ε j + ε u)).card

theorem shiftedBelow_comp {C : Finset ι} {σ : Equiv.Perm ι} (hσ : σ ∈ permsOn C) (s ε : ι → ℝ)
    (u : ι) : shiftedBelow C (s ∘ σ) (ε ∘ σ) u = shiftedBelow C s ε (σ u) := by
  apply card_equiv σ
  intro j
  simp only [mem_filter, Function.comp_apply, mem_iff_of_mem_permsOn hσ]

theorem shiftedAtMost_comp {C : Finset ι} {σ : Equiv.Perm ι} (hσ : σ ∈ permsOn C) (s ε : ι → ℝ)
    (u : ι) : shiftedAtMost C (s ∘ σ) (ε ∘ σ) u = shiftedAtMost C s ε (σ u) := by
  apply card_equiv σ
  intro j
  simp only [mem_filter, Function.comp_apply, mem_iff_of_mem_permsOn hσ]

omit [Fintype ι] [DecidableEq ι] in
theorem shiftedBelow_le_rankLt {C : Finset ι} {s ŝ ε : ι → ℝ}
    (hclose : ∀ j ∈ C, |ŝ j - s j| ≤ ε j) {u : ι} (hu : u ∈ C) :
    shiftedBelow C s ε u ≤ rankLt C ŝ u := by
  refine card_le_card fun j hj => ?_
  rw [mem_filter] at hj ⊢
  have hj' := abs_le.1 (hclose j hj.1)
  have hu' := abs_le.1 (hclose u hu)
  exact ⟨hj.1, by linarith [hj.2]⟩

omit [Fintype ι] [DecidableEq ι] in
theorem rankLe_le_shiftedAtMost {C : Finset ι} {s ŝ ε : ι → ℝ}
    (hclose : ∀ j ∈ C, |ŝ j - s j| ≤ ε j) {u : ι} (hu : u ∈ C) :
    rankLe C ŝ u ≤ shiftedAtMost C s ε u := by
  refine card_le_card fun j hj => ?_
  rw [mem_filter] at hj ⊢
  have hj' := abs_le.1 (hclose j hj.1)
  have hu' := abs_le.1 (hclose u hu)
  exact ⟨hj.1, by linarith [hj.2]⟩

/-- **The PIT from above.** With bounds carried by the units, `P(PIT ≤ k/|C|)` under the estimates
is at most `#{u : shiftedBelow u < k}/|C|`. -/
theorem pit_le_upper {Z : Type*} {C : Finset ι} {t : ι} (ht : t ∈ C) {S : (ι → Z) → ι → ℝ}
    (hS : Equivariant C S) (z : ι → Z) (ŝ : Equiv.Perm ι → ι → ℝ) (ε : ι → ℝ)
    (hclose : ∀ σ ∈ permsOn C, ∀ j ∈ C, |ŝ σ j - S (z ∘ σ) j| ≤ ε (σ j)) (k : ℕ) :
    (∑ σ ∈ permsOn C, pitRankLe C (ŝ σ) t k) / (permsOn C).card ≤
      ((C.filter fun u => shiftedBelow C (S z) ε u < k).card : ℝ) / C.card := by
  classical
  rw [← ratio_permsOn ht (fun u => shiftedBelow C (S z) ε u < k), ← sum_boole]
  have hg : (0 : ℝ) < (permsOn C).card := by exact_mod_cast card_permsOn_pos C
  refine div_le_div_of_nonneg_right (sum_le_sum fun σ hσ => ?_) hg.le
  refine (pitRankLe_le C (ŝ σ) ht k).trans ?_
  have hfun : S (z ∘ σ) = S z ∘ σ := funext fun i => hS σ hσ z i
  have hb : shiftedBelow C (S z) ε (σ t) ≤ rankLt C (ŝ σ) t := by
    rw [← shiftedBelow_comp hσ, ← hfun]
    exact shiftedBelow_le_rankLt (ε := fun j => ε (σ j)) (hclose σ hσ) ht
  split_ifs with h₁ h₂ h₂ <;> first | exact le_rfl | (exfalso; omega) | norm_num

/-- **The PIT from below.** `P(PIT ≤ k/|C|)` under the estimates is at least
`#{u : shiftedAtMost u ≤ k}/|C|`. -/
theorem pit_ge_lower {Z : Type*} {C : Finset ι} {t : ι} (ht : t ∈ C) {S : (ι → Z) → ι → ℝ}
    (hS : Equivariant C S) (z : ι → Z) (ŝ : Equiv.Perm ι → ι → ℝ) (ε : ι → ℝ)
    (hclose : ∀ σ ∈ permsOn C, ∀ j ∈ C, |ŝ σ j - S (z ∘ σ) j| ≤ ε (σ j)) (k : ℕ) :
    ((C.filter fun u => shiftedAtMost C (S z) ε u ≤ k).card : ℝ) / C.card ≤
      (∑ σ ∈ permsOn C, pitRankLe C (ŝ σ) t k) / (permsOn C).card := by
  classical
  rw [← ratio_permsOn ht (fun u => shiftedAtMost C (S z) ε u ≤ k), ← sum_boole]
  have hg : (0 : ℝ) < (permsOn C).card := by exact_mod_cast card_permsOn_pos C
  refine div_le_div_of_nonneg_right (sum_le_sum fun σ hσ => ?_) hg.le
  refine le_trans ?_ (le_pitRankLe C (ŝ σ) ht k)
  have hfun : S (z ∘ σ) = S z ∘ σ := funext fun i => hS σ hσ z i
  have hb : rankLe C (ŝ σ) t ≤ shiftedAtMost C (S z) ε (σ t) := by
    rw [← shiftedAtMost_comp hσ, ← hfun]
    exact rankLe_le_shiftedAtMost (ε := fun j => ε (σ j)) (hclose σ hσ) ht
  split_ifs with h₁ h₂ h₂ <;> first | exact le_rfl | (exfalso; omega) | norm_num

omit [Fintype ι] [DecidableEq ι] in
/-- The two counts bracket `k`: with bounds `≥ 0`, `#{shiftedAtMost ≤ k} ≤ k ≤ #{shiftedBelow < k}`
for `k ≤ |C|`, so the PIT's departure from `k/|C|` is at most their excess over `k`. -/
theorem card_shiftedAtMost_le_le_card_shiftedBelow (C : Finset ι) (s ε : ι → ℝ)
    (hε : ∀ j, 0 ≤ ε j) {k : ℕ} (hk : k ≤ C.card) :
    (C.filter fun u => shiftedAtMost C s ε u ≤ k).card ≤ k ∧
      k ≤ (C.filter fun u => shiftedBelow C s ε u < k).card := by
  constructor
  · refine le_trans (card_le_card fun u hu => ?_) (card_rankLe_le C s k)
    rw [mem_filter] at hu ⊢
    refine ⟨hu.1, le_trans (card_le_card fun j hj => ?_) hu.2⟩
    rw [mem_filter] at hj ⊢
    exact ⟨hj.1, by linarith [hj.2, hε j, hε u]⟩
  · refine le_trans (le_card_rankLt_lt C s hk) (card_le_card fun u hu => ?_)
    rw [mem_filter] at hu ⊢
    refine ⟨hu.1, lt_of_le_of_lt (card_le_card fun j hj => ?_) hu.2⟩
    rw [mem_filter] at hj ⊢
    exact ⟨hj.1, by linarith [hj.2, hε j, hε u]⟩

/-- **The engine's reading for the PIT.** For every ordering, the counts read from its estimates
with doubled shifts bound the exact ones: `#{shiftedBelow < k}` from above and
`#{shiftedAtMost ≤ k}` from below. -/
theorem card_pit_counts_estimated {Z : Type*} {C : Finset ι} {S : (ι → Z) → ι → ℝ}
    (hS : Equivariant C S) (z : ι → Z) {σ : Equiv.Perm ι} (hσ : σ ∈ permsOn C) (ŝ ε : ι → ℝ)
    (hclose : ∀ j ∈ C, |ŝ j - S (z ∘ σ) j| ≤ ε (σ j)) (k : ℕ) :
    (C.filter fun u => shiftedBelow C (S z) ε u < k).card ≤
        (C.filter fun u => shiftedBelow C ŝ (fun j => 2 * ε (σ j)) u < k).card ∧
      (C.filter fun u => shiftedAtMost C ŝ (fun j => 2 * ε (σ j)) u ≤ k).card ≤
        (C.filter fun u => shiftedAtMost C (S z) ε u ≤ k).card := by
  have hfun : S (z ∘ σ) = S z ∘ σ := funext fun i => hS σ hσ z i
  have hB : (C.filter fun u => shiftedBelow C (S z) ε u < k).card =
      (C.filter fun u => shiftedBelow C (S (z ∘ σ)) (ε ∘ σ) u < k).card := by
    symm
    apply card_equiv σ
    intro u
    simp only [mem_filter, hfun, shiftedBelow_comp hσ, mem_iff_of_mem_permsOn hσ]
  have hA : (C.filter fun u => shiftedAtMost C (S z) ε u ≤ k).card =
      (C.filter fun u => shiftedAtMost C (S (z ∘ σ)) (ε ∘ σ) u ≤ k).card := by
    symm
    apply card_equiv σ
    intro u
    simp only [mem_filter, hfun, shiftedAtMost_comp hσ, mem_iff_of_mem_permsOn hσ]
  rw [hB, hA]
  constructor
  · refine card_le_card fun u hu => ?_
    rw [mem_filter] at hu ⊢
    refine ⟨hu.1, lt_of_le_of_lt (card_le_card fun j hj => ?_) hu.2⟩
    rw [mem_filter] at hj ⊢
    have hj' := abs_le.1 (hclose j hj.1)
    have hu' := abs_le.1 (hclose u hu.1)
    exact ⟨hj.1, by simp only [Function.comp_apply]; linarith [hj.2]⟩
  · refine card_le_card fun u hu => ?_
    rw [mem_filter] at hu ⊢
    refine ⟨hu.1, le_trans (card_le_card fun j hj => ?_) hu.2⟩
    rw [mem_filter] at hj ⊢
    have hj' := abs_le.1 (hclose j hj.1)
    have hu' := abs_le.1 (hclose u hu.1)
    exact ⟨hj.1, by simp only [Function.comp_apply] at hj ⊢; linarith [hj.2]⟩

/-! #### The carried bound that leaves out the degenerate turns -/

/-- **An average under a chance of failure, from above.** The draw `w` reads `f w σ ≤ 1` on the
ordering `σ`. On an ordering whose test unit is kept (not `out w`) and whose bounds hold
(`hold w σ`), the reading is at most `0` unless the test unit lands on a `good` unit. The average
reading is at most the share of units left out or good, plus `ν`, the weighted share of kept
orderings whose bounds fail. -/
theorem avg_le_of_fail {W : Type*} [Fintype W] {C : Finset ι} {t : ι} (ht : t ∈ C)
    (p : W → ℝ) (hp0 : ∀ w, 0 ≤ p w) (out good : W → ι → Prop) [∀ w, DecidablePred (out w)]
    [∀ w, DecidablePred (good w)] (hold : W → Equiv.Perm ι → Prop)
    [∀ w, DecidablePred (hold w)] (f : W → Equiv.Perm ι → ℝ)
    (hf1 : ∀ w, ∀ σ ∈ permsOn C, f w σ ≤ 1)
    (h : ∀ w, ∀ σ ∈ permsOn C, ¬ out w (σ t) → hold w σ → good w (σ t) ∨ f w σ ≤ 0) {ν : ℝ}
    (hfail : ∑ w, p w * (((permsOn C).filter fun σ : Equiv.Perm ι =>
        ¬ out w (σ t) ∧ ¬ hold w σ).card : ℝ) ≤ ν * (permsOn C).card) :
    ∑ w, p w * ((∑ σ ∈ permsOn C, f w σ) / (permsOn C).card) ≤
      ∑ w, p w * (((C.filter fun u => out w u ∨ good w u).card : ℝ) / C.card) + ν := by
  have hg : (0 : ℝ) < (permsOn C).card := by exact_mod_cast card_permsOn_pos C
  have hw : ∀ w, (∑ σ ∈ permsOn C, f w σ) / (permsOn C).card ≤
      ((C.filter fun u => out w u ∨ good w u).card : ℝ) / C.card +
        (((permsOn C).filter fun σ : Equiv.Perm ι => ¬ out w (σ t) ∧ ¬ hold w σ).card : ℝ) /
          (permsOn C).card := by
    intro w
    rw [← ratio_permsOn ht (fun u => out w u ∨ good w u), ← add_div, ← sum_boole, ← sum_boole,
      ← sum_add_distrib]
    refine div_le_div_of_nonneg_right (sum_le_sum fun σ hσ => ?_) hg.le
    have h1 := hf1 w σ hσ
    split_ifs with hA hB hB
    · linarith
    · linarith
    · linarith
    · simp only [not_or, not_and, not_not] at hA hB
      rcases h w σ hσ hA.1 (hB hA.1) with hgd | hf
      · exact absurd hgd hA.2
      · linarith
  have hsum := sum_le_sum fun w (_ : w ∈ univ) => mul_le_mul_of_nonneg_left (hw w) (hp0 w)
  simp only [mul_add, sum_add_distrib] at hsum
  have hν : ∑ w, p w * ((((permsOn C).filter fun σ : Equiv.Perm ι =>
      ¬ out w (σ t) ∧ ¬ hold w σ).card : ℝ) / (permsOn C).card) ≤ ν :=
    weighted_share_le C p (fun w σ => ¬ out w (σ t) ∧ ¬ hold w σ) hfail
  linarith

/-- **The same from below.** A reading `f w σ ≥ 0` that is at least `1` on a kept ordering whose
bounds hold and whose test unit lands on a good unit averages at least the share of kept good
units, less `ν`. -/
theorem avg_ge_of_fail {W : Type*} [Fintype W] {C : Finset ι} {t : ι} (ht : t ∈ C)
    (p : W → ℝ) (hp0 : ∀ w, 0 ≤ p w) (out good : W → ι → Prop) [∀ w, DecidablePred (out w)]
    [∀ w, DecidablePred (good w)] (hold : W → Equiv.Perm ι → Prop)
    [∀ w, DecidablePred (hold w)] (f : W → Equiv.Perm ι → ℝ)
    (hf0 : ∀ w, ∀ σ ∈ permsOn C, 0 ≤ f w σ)
    (h : ∀ w, ∀ σ ∈ permsOn C, ¬ out w (σ t) → hold w σ → good w (σ t) → 1 ≤ f w σ) {ν : ℝ}
    (hfail : ∑ w, p w * (((permsOn C).filter fun σ : Equiv.Perm ι =>
        ¬ out w (σ t) ∧ ¬ hold w σ).card : ℝ) ≤ ν * (permsOn C).card) :
    ∑ w, p w * (((C.filter fun u => ¬ out w u ∧ good w u).card : ℝ) / C.card) - ν ≤
      ∑ w, p w * ((∑ σ ∈ permsOn C, f w σ) / (permsOn C).card) := by
  have hg : (0 : ℝ) < (permsOn C).card := by exact_mod_cast card_permsOn_pos C
  have hw : ∀ w, ((C.filter fun u => ¬ out w u ∧ good w u).card : ℝ) / C.card ≤
      (∑ σ ∈ permsOn C, f w σ) / (permsOn C).card +
        (((permsOn C).filter fun σ : Equiv.Perm ι => ¬ out w (σ t) ∧ ¬ hold w σ).card : ℝ) /
          (permsOn C).card := by
    intro w
    rw [← ratio_permsOn ht (fun u => ¬ out w u ∧ good w u), ← add_div, ← sum_boole,
      ← sum_boole, ← sum_add_distrib]
    refine div_le_div_of_nonneg_right (sum_le_sum fun σ hσ => ?_) hg.le
    have h0 := hf0 w σ hσ
    split_ifs with hA hB hB
    · linarith
    · simp only [not_and, not_not] at hB
      linarith [h w σ hσ hA.1 (hB hA.1) hA.2]
    · linarith
    · linarith
  have hsum := sum_le_sum fun w (_ : w ∈ univ) => mul_le_mul_of_nonneg_left (hw w) (hp0 w)
  simp only [mul_add, sum_add_distrib] at hsum
  have hν : ∑ w, p w * ((((permsOn C).filter fun σ : Equiv.Perm ι =>
      ¬ out w (σ t) ∧ ¬ hold w σ).card : ℝ) / (permsOn C).card) ≤ ν :=
    weighted_share_le C p (fun w σ => ¬ out w (σ t) ∧ ¬ hold w σ) hfail
  linarith

/-- **The margin, leaving out the degenerate turns.** `W` carries the draws, with weights `p`. For
the draw `w`, `D w` is the set of turns whose weights degenerate (read once for the data set, by
`k̂` or the ESS, so the same for every ordering), `ε w` the bounds carried over the other turns
only, and `ν` bounds the weighted share of kept orderings on which some unit's bound fails. Then
coverage is at least `1 − |D ∪ R ∪ M_D|/|C| − ν`, with `R` the exact rule's refusals and `M_D` its
fragile units under the kept bounds. Since `|R| ≤ k`, this is at least
`1 − (k + |M_D| + |D|)/|C| − ν`. -/
theorem coverage_close_excluding {W Z : Type*} [Fintype W] {C : Finset ι} {t : ι} (ht : t ∈ C)
    {S : (ι → Z) → ι → ℝ} (hS : Equivariant C S) (z : ι → Z) (p : W → ℝ) (hp0 : ∀ w, 0 ≤ p w)
    (hp1 : ∑ w, p w = 1) (ŝ : W → Equiv.Perm ι → ι → ℝ) (ε : W → ι → ℝ) (D : W → Finset ι)
    {ν : ℝ} (hfail : ∑ w, p w * (((permsOn C).filter fun σ : Equiv.Perm ι => σ t ∉ D w ∧
        ¬ ∀ j ∈ C, |ŝ w σ j - S (z ∘ σ) j| ≤ ε w (σ j)).card : ℝ) ≤ ν * (permsOn C).card)
    (k : ℕ) :
    1 - ∑ w, p w * (((C.filter fun u => u ∈ D w ∨ ¬ Accepts k C (S z) u ∨
        Fragile k C (S z) (ε w) u).card : ℝ) / C.card) - ν ≤
      ∑ w, p w * ((((permsOn C).filter fun σ : Equiv.Perm ι =>
        Accepts k C (ŝ w σ) t).card : ℝ) / (permsOn C).card) := by
  have key := coverage_ge_of_good_fail ht p hp0
    (fun w u => ¬ (u ∈ D w ∨ ¬ Accepts k C (S z) u ∨ Fragile k C (S z) (ε w) u))
    (fun w σ => ∀ j ∈ C, |ŝ w σ j - S (z ∘ σ) j| ≤ ε w (σ j))
    (fun w σ => Accepts k C (ŝ w σ) t) (fun w σ hσ hg hh => by
      simp only [not_or, not_not] at hg
      have hfun : S (z ∘ σ) = S z ∘ σ := funext fun i => hS σ hσ z i
      have hacc : Accepts k C (S (z ∘ σ)) t := by
        rw [hfun]
        simp only [Accepts, rankLe_comp hσ]
        exact hg.2.1
      rcases accepts_or_fragile (ε := fun j => ε w (σ j)) hh ht hacc with h | h
      · exact h
      · exfalso
        apply hg.2.2
        rw [hfun] at h
        exact (fragile_comp hσ (S z) (ε w) t).1 h)
    (weighted_fail_mono C p hp0 _ _ (by
      intro w σ h
      exact ⟨fun hD => h.1 (Or.inl hD), h.2⟩) hfail)
  rw [weighted_share_not ⟨t, ht⟩ p hp1] at key
  exact key

/-- **The inflated rule, leaving out the degenerate turns.** With bounds by ordering and position,
trusted on the kept orderings with the same `ν`, coverage is at least `1 − |D ∪ R|/|C| − ν`. -/
theorem coverage_inflated_excluding {W Z : Type*} [Fintype W] {C : Finset ι} {t : ι}
    (ht : t ∈ C) {S : (ι → Z) → ι → ℝ} (hS : Equivariant C S) (z : ι → Z) (p : W → ℝ)
    (hp0 : ∀ w, 0 ≤ p w) (hp1 : ∑ w, p w = 1) (ŝ ε : W → Equiv.Perm ι → ι → ℝ)
    (D : W → Finset ι) {ν : ℝ}
    (hfail : ∑ w, p w * (((permsOn C).filter fun σ : Equiv.Perm ι => σ t ∉ D w ∧
        ¬ ∀ j ∈ C, |ŝ w σ j - S (z ∘ σ) j| ≤ ε w σ j).card : ℝ) ≤ ν * (permsOn C).card)
    (k : ℕ) :
    1 - ∑ w, p w * (((C.filter fun u => u ∈ D w ∨ ¬ Accepts k C (S z) u).card : ℝ) / C.card) -
        ν ≤
      ∑ w, p w * ((((permsOn C).filter fun σ : Equiv.Perm ι =>
        AcceptsInflated k C (ŝ w σ) (ε w σ) t).card : ℝ) / (permsOn C).card) := by
  have key := coverage_ge_of_good_fail ht p hp0 (fun w u => ¬ (u ∈ D w ∨ ¬ Accepts k C (S z) u))
    (fun w σ => ∀ j ∈ C, |ŝ w σ j - S (z ∘ σ) j| ≤ ε w σ j)
    (fun w σ => AcceptsInflated k C (ŝ w σ) (ε w σ) t) (fun w σ hσ hg hh => by
      simp only [not_or, not_not] at hg
      have hfun : S (z ∘ σ) = S z ∘ σ := funext fun i => hS σ hσ z i
      have hacc : Accepts k C (S (z ∘ σ)) t := by
        rw [hfun]
        simp only [Accepts, rankLe_comp hσ]
        exact hg.2
      exact acceptsInflated_of_accepts hh ht hacc)
    (weighted_fail_mono C p hp0 _ _ (by
      intro w σ h
      exact ⟨fun hD => h.1 (Or.inl hD), h.2⟩) hfail)
  rw [weighted_share_not ⟨t, ht⟩ p hp1] at key
  exact key

/-- **The equi-tailed set, leaving out the degenerate turns.** Coverage is at least
`1 − |D ∪ R ∪ M_lo ∪ M_hi|/|C| − ν`, with `R` the exact set's refusals (at most `2k`) and the
tails' bands read under the kept bounds. -/
theorem coverage_close_twoSided_excluding {W Z : Type*} [Fintype W] {C : Finset ι} {t : ι}
    (ht : t ∈ C) {S : (ι → Z) → ι → ℝ} (hS : Equivariant C S) (z : ι → Z) (p : W → ℝ)
    (hp0 : ∀ w, 0 ≤ p w) (hp1 : ∑ w, p w = 1) (ŝ : W → Equiv.Perm ι → ι → ℝ) (ε : W → ι → ℝ)
    (D : W → Finset ι) {ν : ℝ}
    (hfail : ∑ w, p w * (((permsOn C).filter fun σ : Equiv.Perm ι => σ t ∉ D w ∧
        ¬ ∀ j ∈ C, |ŝ w σ j - S (z ∘ σ) j| ≤ ε w (σ j)).card : ℝ) ≤ ν * (permsOn C).card)
    (k : ℕ) :
    1 - ∑ w, p w * (((C.filter fun u => u ∈ D w ∨ ¬ AcceptsTwoSided k C (S z) u ∨
        FragileTwoSided k C (S z) (ε w) u).card : ℝ) / C.card) - ν ≤
      ∑ w, p w * ((((permsOn C).filter fun σ : Equiv.Perm ι =>
        AcceptsTwoSided k C (ŝ w σ) t).card : ℝ) / (permsOn C).card) := by
  have key := coverage_ge_of_good_fail ht p hp0
    (fun w u => ¬ (u ∈ D w ∨ ¬ AcceptsTwoSided k C (S z) u ∨ FragileTwoSided k C (S z) (ε w) u))
    (fun w σ => ∀ j ∈ C, |ŝ w σ j - S (z ∘ σ) j| ≤ ε w (σ j))
    (fun w σ => AcceptsTwoSided k C (ŝ w σ) t) (fun w σ hσ hg hh => by
      simp only [not_or, not_not] at hg
      have hfun : S (z ∘ σ) = S z ∘ σ := funext fun i => hS σ hσ z i
      have hacc : AcceptsTwoSided k C (S (z ∘ σ)) t := by
        rw [hfun, acceptsTwoSided_comp hσ]
        exact hg.2.1
      rcases acceptsTwoSided_or_fragile (ε := fun j => ε w (σ j)) hh ht hacc with h | h
      · exact h
      · exfalso
        apply hg.2.2
        rw [hfun] at h
        exact (fragileTwoSided_comp hσ (S z) (ε w) t).1 h)
    (weighted_fail_mono C p hp0 _ _ (by
      intro w σ h
      exact ⟨fun hD => h.1 (Or.inl hD), h.2⟩) hfail)
  rw [weighted_share_not ⟨t, ht⟩ p hp1] at key
  exact key

/-- **The inflated equi-tailed set, leaving out the degenerate turns:** at least
`1 − |D ∪ R|/|C| − ν`, so at least `1 − (2k + |D|)/|C| − ν`. -/
theorem coverage_inflated_twoSided_excluding {W Z : Type*} [Fintype W] {C : Finset ι} {t : ι}
    (ht : t ∈ C) {S : (ι → Z) → ι → ℝ} (hS : Equivariant C S) (z : ι → Z) (p : W → ℝ)
    (hp0 : ∀ w, 0 ≤ p w) (hp1 : ∑ w, p w = 1) (ŝ ε : W → Equiv.Perm ι → ι → ℝ)
    (D : W → Finset ι) {ν : ℝ}
    (hfail : ∑ w, p w * (((permsOn C).filter fun σ : Equiv.Perm ι => σ t ∉ D w ∧
        ¬ ∀ j ∈ C, |ŝ w σ j - S (z ∘ σ) j| ≤ ε w σ j).card : ℝ) ≤ ν * (permsOn C).card)
    (k : ℕ) :
    1 - ∑ w, p w * (((C.filter fun u => u ∈ D w ∨ ¬ AcceptsTwoSided k C (S z) u).card : ℝ) /
        C.card) - ν ≤
      ∑ w, p w * ((((permsOn C).filter fun σ : Equiv.Perm ι =>
        AcceptsInflated k C (ŝ w σ) (ε w σ) t ∧
          AcceptsInflated k C (fun j => -ŝ w σ j) (ε w σ) t).card : ℝ) / (permsOn C).card) := by
  have key := coverage_ge_of_good_fail ht p hp0
    (fun w u => ¬ (u ∈ D w ∨ ¬ AcceptsTwoSided k C (S z) u))
    (fun w σ => ∀ j ∈ C, |ŝ w σ j - S (z ∘ σ) j| ≤ ε w σ j)
    (fun w σ => AcceptsInflated k C (ŝ w σ) (ε w σ) t ∧
      AcceptsInflated k C (fun j => -ŝ w σ j) (ε w σ) t) (fun w σ hσ hg hh => by
      simp only [not_or, not_not] at hg
      have hfun : S (z ∘ σ) = S z ∘ σ := funext fun i => hS σ hσ z i
      have hacc : AcceptsTwoSided k C (S (z ∘ σ)) t := by
        rw [hfun, acceptsTwoSided_comp hσ]
        exact hg.2
      exact ⟨acceptsInflated_of_accepts hh ht hacc.1,
        acceptsInflated_of_accepts (close_neg hh) ht hacc.2⟩)
    (weighted_fail_mono C p hp0 _ _ (by
      intro w σ h
      exact ⟨fun hD => h.1 (Or.inl hD), h.2⟩) hfail)
  rw [weighted_share_not ⟨t, ht⟩ p hp1] at key
  exact key

/-- **The PIT, leaving out the degenerate turns.** With bounds carried over the kept turns, failing
on a weighted share `ν` of kept orderings, `P(PIT ≤ k/|C|)` lies between
`#{u ∉ D : shiftedAtMost u ≤ k}/|C| − ν` and `#{u ∈ D or shiftedBelow u < k}/|C| + ν`. -/
theorem pit_excluding {W Z : Type*} [Fintype W] {C : Finset ι} {t : ι} (ht : t ∈ C)
    {S : (ι → Z) → ι → ℝ} (hS : Equivariant C S) (z : ι → Z) (p : W → ℝ) (hp0 : ∀ w, 0 ≤ p w)
    (ŝ : W → Equiv.Perm ι → ι → ℝ) (ε : W → ι → ℝ) (D : W → Finset ι) {ν : ℝ}
    (hfail : ∑ w, p w * (((permsOn C).filter fun σ : Equiv.Perm ι => σ t ∉ D w ∧
        ¬ ∀ j ∈ C, |ŝ w σ j - S (z ∘ σ) j| ≤ ε w (σ j)).card : ℝ) ≤ ν * (permsOn C).card)
    (k : ℕ) :
    ∑ w, p w * (((C.filter fun u => u ∉ D w ∧ shiftedAtMost C (S z) (ε w) u ≤ k).card : ℝ) /
        C.card) - ν ≤
      ∑ w, p w * ((∑ σ ∈ permsOn C, pitRankLe C (ŝ w σ) t k) / (permsOn C).card) ∧
    ∑ w, p w * ((∑ σ ∈ permsOn C, pitRankLe C (ŝ w σ) t k) / (permsOn C).card) ≤
      ∑ w, p w * (((C.filter fun u => u ∈ D w ∨ shiftedBelow C (S z) (ε w) u < k).card : ℝ) /
        C.card) + ν := by
  constructor
  · refine avg_ge_of_fail ht p hp0 (fun w u => u ∈ D w)
      (fun w u => shiftedAtMost C (S z) (ε w) u ≤ k)
      (fun w σ => ∀ j ∈ C, |ŝ w σ j - S (z ∘ σ) j| ≤ ε w (σ j))
      (fun w σ => pitRankLe C (ŝ w σ) t k) (fun w σ _ => pitRankLe_nonneg _ _ _ _)
      (fun w σ hσ _ hh hg => ?_) hfail
    have hfun : S (z ∘ σ) = S z ∘ σ := funext fun i => hS σ hσ z i
    have hb : rankLe C (ŝ w σ) t ≤ shiftedAtMost C (S z) (ε w) (σ t) := by
      rw [← shiftedAtMost_comp hσ, ← hfun]
      exact rankLe_le_shiftedAtMost (ε := fun j => ε w (σ j)) hh ht
    have hle := le_pitRankLe C (ŝ w σ) ht k
    rw [if_pos (hb.trans hg)] at hle
    exact hle
  · refine avg_le_of_fail ht p hp0 (fun w u => u ∈ D w)
      (fun w u => shiftedBelow C (S z) (ε w) u < k)
      (fun w σ => ∀ j ∈ C, |ŝ w σ j - S (z ∘ σ) j| ≤ ε w (σ j))
      (fun w σ => pitRankLe C (ŝ w σ) t k) (fun w σ _ => ?_) (fun w σ hσ _ hh => ?_) hfail
    · exact (pitRankLe_le C (ŝ w σ) ht k).trans (by split_ifs <;> norm_num)
    · by_cases hg : shiftedBelow C (S z) (ε w) (σ t) < k
      · exact Or.inl hg
      · refine Or.inr ?_
        have hfun : S (z ∘ σ) = S z ∘ σ := funext fun i => hS σ hσ z i
        have hb : shiftedBelow C (S z) (ε w) (σ t) ≤ rankLt C (ŝ w σ) t := by
          rw [← shiftedBelow_comp hσ, ← hfun]
          exact shiftedBelow_le_rankLt (ε := fun j => ε w (σ j)) hh ht
        have hlt : ¬ rankLt C (ŝ w σ) t < k := by omega
        have hle := pitRankLe_le C (ŝ w σ) ht k
        rw [if_neg hlt] at hle
        exact hle

/-! #### What the engine reads, leaving out the degenerate turns -/

omit [Fintype ι] in
/-- **The share, read from estimates.** From estimates within the kept bounds on one ordering
(units at their positions), the share in `coverage_close_excluding` is at most `k`, plus the band
read from the estimates, plus the units of `D` that the inflated rule accepts. A degenerate turn
that the inflated rule refuses is refused by the exact rule too, so it costs nothing. -/
theorem card_excluded_le_estimated {k : ℕ} {C : Finset ι} {s ŝ ε : ι → ℝ}
    (hclose : ∀ j ∈ C, |ŝ j - s j| ≤ ε j) (D : Finset ι) :
    (C.filter fun u => u ∈ D ∨ ¬ Accepts k C s u ∨ Fragile k C s ε u).card ≤
      k + (C.filter fun u => FragileEstimated k C ŝ ε u).card +
        (C.filter fun u => u ∈ D ∧ AcceptsInflated k C ŝ ε u).card := by
  have hsub : (C.filter fun u => u ∈ D ∨ ¬ Accepts k C s u ∨ Fragile k C s ε u) ⊆
      ((C.filter fun u => rankLe C s u ≤ k) ∪ C.filter fun u => FragileEstimated k C ŝ ε u) ∪
        C.filter fun u => u ∈ D ∧ AcceptsInflated k C ŝ ε u := by
    intro u hu
    rw [mem_filter] at hu
    simp only [mem_union, mem_filter]
    by_cases hacc : Accepts k C s u
    · rcases hu.2 with hD | hna | hfr
      · exact Or.inr ⟨hu.1, hD, acceptsInflated_of_accepts hclose hu.1 hacc⟩
      · exact absurd hacc hna
      · exact Or.inl (Or.inr ⟨hu.1, fragileEstimated_of_fragile hclose hu.1 hfr⟩)
    · exact Or.inl (Or.inl ⟨hu.1, not_lt.1 hacc⟩)
  have h1 := card_rankLe_le C s k
  have h2 := (card_le_card hsub).trans ((card_union_le _ _).trans
    (Nat.add_le_add_right (card_union_le _ _) _))
  omega

omit [Fintype ι] in
/-- **The same for the equi-tailed set:** at most `2k`, plus the union band read from the
estimates, plus the units of `D` the inflated set accepts. -/
theorem card_excluded_le_estimated_twoSided {k : ℕ} {C : Finset ι} {s ŝ ε : ι → ℝ}
    (hclose : ∀ j ∈ C, |ŝ j - s j| ≤ ε j) (D : Finset ι) :
    (C.filter fun u => u ∈ D ∨ ¬ AcceptsTwoSided k C s u ∨ FragileTwoSided k C s ε u).card ≤
      2 * k + (C.filter fun u => FragileEstimated k C ŝ ε u ∨
          FragileEstimated k C (fun j => -ŝ j) ε u).card +
        (C.filter fun u => u ∈ D ∧ AcceptsInflated k C ŝ ε u ∧
          AcceptsInflated k C (fun j => -ŝ j) ε u).card := by
  have hclose' := close_neg hclose
  have hsub : (C.filter fun u => u ∈ D ∨ ¬ AcceptsTwoSided k C s u ∨
      FragileTwoSided k C s ε u) ⊆
      ((C.filter fun u => ¬ AcceptsTwoSided k C s u) ∪ C.filter fun u =>
        FragileEstimated k C ŝ ε u ∨ FragileEstimated k C (fun j => -ŝ j) ε u) ∪
        C.filter fun u => u ∈ D ∧ AcceptsInflated k C ŝ ε u ∧
          AcceptsInflated k C (fun j => -ŝ j) ε u := by
    intro u hu
    rw [mem_filter] at hu
    simp only [mem_union, mem_filter]
    by_cases hacc : AcceptsTwoSided k C s u
    · rcases hu.2 with hD | hna | hfr
      · exact Or.inr ⟨hu.1, hD, acceptsInflated_of_accepts hclose hu.1 hacc.1,
          acceptsInflated_of_accepts hclose' hu.1 hacc.2⟩
      · exact absurd hacc hna
      · refine Or.inl (Or.inr ⟨hu.1, ?_⟩)
        rcases hfr with h | h
        · exact Or.inl (fragileEstimated_of_fragile hclose hu.1 h)
        · exact Or.inr (fragileEstimated_of_fragile hclose' hu.1 h)
    · exact Or.inl (Or.inl ⟨hu.1, hacc⟩)
  have h1 := card_not_acceptsTwoSided_le C s k
  have h2 := (card_le_card hsub).trans ((card_union_le _ _).trans
    (Nat.add_le_add_right (card_union_le _ _) _))
  omega

omit [Fintype ι] [DecidableEq ι] in
/-- Read from estimates with a shift at least doubled, fewer units sit below `u`. -/
theorem shiftedBelow_estimated_le {C : Finset ι} {s ŝ ε δ : ι → ℝ}
    (hclose : ∀ j ∈ C, |ŝ j - s j| ≤ ε j) (hδ : ∀ j ∈ C, 2 * ε j ≤ δ j) {u : ι} (hu : u ∈ C) :
    shiftedBelow C ŝ δ u ≤ shiftedBelow C s ε u := by
  refine card_le_card fun j hj => ?_
  rw [mem_filter] at hj ⊢
  have hj' := abs_le.1 (hclose j hj.1)
  have hu' := abs_le.1 (hclose u hu)
  exact ⟨hj.1, by linarith [hj.2, hδ j hj.1, hδ u hu]⟩

omit [Fintype ι] [DecidableEq ι] in
/-- Read from estimates with a shift at least doubled, more units sit at most above `u`. -/
theorem shiftedAtMost_le_estimated {C : Finset ι} {s ŝ ε δ : ι → ℝ}
    (hclose : ∀ j ∈ C, |ŝ j - s j| ≤ ε j) (hδ : ∀ j ∈ C, 2 * ε j ≤ δ j) {u : ι} (hu : u ∈ C) :
    shiftedAtMost C s ε u ≤ shiftedAtMost C ŝ δ u := by
  refine card_le_card fun j hj => ?_
  rw [mem_filter] at hj ⊢
  have hj' := abs_le.1 (hclose j hj.1)
  have hu' := abs_le.1 (hclose u hu)
  exact ⟨hj.1, by linarith [hj.2, hδ j hj.1, hδ u hu]⟩

omit [Fintype ι] in
/-- **The PIT's band with `D`, read from estimates.** The counts in `pit_excluding`, read from one
ordering's estimates with the shifts doubled, bound the exact ones from the safe side. -/
theorem card_pit_excluding_estimated {C : Finset ι} {s ŝ ε : ι → ℝ}
    (hclose : ∀ j ∈ C, |ŝ j - s j| ≤ ε j) (D : Finset ι) (k : ℕ) :
    (C.filter fun u => u ∈ D ∨ shiftedBelow C s ε u < k).card ≤
        (C.filter fun u => u ∈ D ∨ shiftedBelow C ŝ (fun j => 2 * ε j) u < k).card ∧
      (C.filter fun u => u ∉ D ∧ shiftedAtMost C ŝ (fun j => 2 * ε j) u ≤ k).card ≤
        (C.filter fun u => u ∉ D ∧ shiftedAtMost C s ε u ≤ k).card := by
  have hδ : ∀ j ∈ C, 2 * ε j ≤ (fun j => 2 * ε j) j := fun j _ => le_rfl
  constructor
  · refine card_le_card fun u hu => ?_
    rw [mem_filter] at hu ⊢
    exact ⟨hu.1, hu.2.imp id fun h =>
      lt_of_le_of_lt (shiftedBelow_estimated_le hclose hδ hu.1) h⟩
  · refine card_le_card fun u hu => ?_
    rw [mem_filter] at hu ⊢
    exact ⟨hu.1, hu.2.1, le_trans (shiftedAtMost_le_estimated hclose hδ hu.1) hu.2.2⟩

end FullStepError

/-! ### Slice 12: Hoff's score and the latent posterior -/

section ShiftedBounds

variable {ι : Type*} [DecidableEq ι]

omit [DecidableEq ι] in
/-- Ranks see only differences: one constant added to every score changes no rank. -/
theorem rankLe_add_const (C : Finset ι) (s : ι → ℝ) (c : ℝ) (u : ι) :
    rankLe C (fun j => s j + c) u = rankLe C s u := by
  unfold rankLe
  congr 1
  exact filter_congr fun j _ => add_le_add_iff_right c

omit [DecidableEq ι] in
theorem accepts_add_const_iff (k : ℕ) (C : Finset ι) (s : ι → ℝ) (c : ℝ) (u : ι) :
    Accepts k C (fun j => s j + c) u ↔ Accepts k C s u := by
  unfold Accepts
  rw [rankLe_add_const]

omit [DecidableEq ι] in
theorem acceptsInflated_add_const_iff (k : ℕ) (C : Finset ι) (s ε : ι → ℝ) (c : ℝ) (u : ι) :
    AcceptsInflated k C (fun j => s j + c) ε u ↔ AcceptsInflated k C s ε u := by
  unfold AcceptsInflated
  rw [filter_congr fun j _ => show s j + c - ε j ≤ s u + c + ε u ↔ s j - ε j ≤ s u + ε u by
    constructor <;> intro h <;> linarith]

variable [Fintype ι]

/-- **CB-3 with bounds up to a common shift per ordering.** The estimates need only lie within
their bounds of the exact scores after one constant `c σ` is added to every unit's estimate in the
ordering `σ`. A factor common to every unit of an ordering, such as a shared normalizer on the log
scale, so drops out of the bounds. The same transfer holds for every form of CB-3 and slice 10,
since the rules see only differences. -/
theorem coverage_close_shift {Z : Type*} {C : Finset ι} {t : ι} (ht : t ∈ C)
    {S : (ι → Z) → ι → ℝ} (hS : Equivariant C S) (z : ι → Z) (ŝ : Equiv.Perm ι → ι → ℝ)
    (c : Equiv.Perm ι → ℝ) (ε : ι → ℝ)
    (hclose : ∀ σ ∈ permsOn C, ∀ j ∈ C, |ŝ σ j + c σ - S (z ∘ σ) j| ≤ ε (σ j)) (k : ℕ) :
    1 - ((k : ℝ) + (C.filter fun u => Fragile k C (S z) ε u).card) / C.card ≤
      (((permsOn C).filter fun σ : Equiv.Perm ι => Accepts k C (ŝ σ) t).card : ℝ) /
        (permsOn C).card := by
  have h := coverage_close ht hS z (fun σ j => ŝ σ j + c σ) ε hclose k
  simp only [accepts_add_const_iff] at h
  exact h

/-- **The inflated rule with bounds up to a common shift per ordering.** -/
theorem coverage_inflated_shift {Z : Type*} {C : Finset ι} {t : ι} (ht : t ∈ C)
    {S : (ι → Z) → ι → ℝ} (hS : Equivariant C S) (z : ι → Z) (ŝ ε : Equiv.Perm ι → ι → ℝ)
    (c : Equiv.Perm ι → ℝ)
    (hclose : ∀ σ ∈ permsOn C, ∀ j ∈ C, |ŝ σ j + c σ - S (z ∘ σ) j| ≤ ε σ j) (k : ℕ) :
    1 - (k : ℝ) / C.card ≤
      (((permsOn C).filter fun σ : Equiv.Perm ι => AcceptsInflated k C (ŝ σ) (ε σ) t).card : ℝ) /
        (permsOn C).card := by
  have h := coverage_inflated ht hS z (fun σ j => ŝ σ j + c σ) ε hclose k
  simp only [acceptsInflated_add_const_iff] at h
  exact h

end ShiftedBounds

section Hoff

variable {Θ V Z : Type*} [Fintype Θ] [Fintype V]
variable (prior : Θ → ℝ) (q : Θ → V → ℝ) (f : Θ → V → Z → ℝ)
variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- The unnormalized posterior of `θ` given every unit but `i`, each latent integrated. -/
noncomputable def postUExcept (z : ι → Z) (i : ι) (θ : Θ) : ℝ :=
  prior θ * ∏ j ∈ univ.erase i, integrated q f θ (z j)

/-- **Hoff's score**, every latent integrated (CB-4(b)'s reading): unit `i`'s predictive density
given every other unit, `c_i = p(z_i | z_{-i})`. Hoff 2023, §4.1. -/
noncomputable def hoffScore (z : ι → Z) (i : ι) : ℝ :=
  (∑ θ, postUExcept prior q f z i θ * integrated q f θ (z i)) / ∑ θ, postUExcept prior q f z i θ

omit [Fintype Θ] in
theorem postUExcept_mul (z : ι → Z) (i : ι) (θ : Θ) :
    postUExcept prior q f z i θ * integrated q f θ (z i) = postU prior q f z θ := by
  rw [postUExcept, postU, ← mul_prod_erase univ (fun j => integrated q f θ (z j)) (mem_univ i)]
  ring

/-- `c_i = p(z) / p(z_{-i})`: the numerator is the same for every unit. -/
theorem hoffScore_eq (z : ι → Z) (i : ι) :
    hoffScore prior q f z i = (∑ θ, postU prior q f z θ) / ∑ θ, postUExcept prior q f z i θ := by
  simp only [hoffScore, postUExcept_mul]

omit [Fintype Θ] in
theorem postUExcept_comp (z : ι → Z) (σ : Equiv.Perm ι) (i : ι) (θ : Θ) :
    postUExcept prior q f (z ∘ σ) i θ = postUExcept prior q f z (σ i) θ := by
  unfold postUExcept
  congr 1
  exact prod_equiv σ (fun j => by simp) fun j _ => rfl

theorem hoffScore_comp (z : ι → Z) (σ : Equiv.Perm ι) (i : ι) :
    hoffScore prior q f (z ∘ σ) i = hoffScore prior q f z (σ i) := by
  simp only [hoffScore_eq, postU_comp, postUExcept_comp]

/-- **Hoff's score is equivariant**, so CB-1 covers with it (`hoff_coverage`), ties counted for
the candidate. Hoff's strict rule does not cover under ties (the CB-4 finding). -/
theorem hoffScore_equivariant (C : Finset ι) : Equivariant C (hoffScore prior q f) :=
  fun σ _ z i => hoffScore_comp prior q f z σ i

theorem hoff_coverage {C : Finset ι} {t : ι} (ht : t ∈ C) (z : ι → Z) (k : ℕ) :
    1 - (k : ℝ) / C.card ≤
      (((permsOn C).filter fun σ : Equiv.Perm ι =>
        Accepts k C (hoffScore prior q f (z ∘ σ)) t).card : ℝ) / (permsOn C).card :=
  coverage_orderings ht (hoffScore_equivariant prior q f C) z k

/-! #### Hoff's score from the `n`-unit draws -/

omit [Fintype Θ] in
theorem postUExcept_none {ι₀ : Type*} [Fintype ι₀] [DecidableEq ι₀] (z : Option ι₀ → Z)
    (θ : Θ) : postUExcept prior q f z none θ = postU prior q f (fun j => z (some j)) θ := by
  have he : (univ : Finset (Option ι₀)).erase none = univ.map Function.Embedding.some := by
    ext o
    cases o <;> simp
  rw [postUExcept, postU, he, prod_map]
  rfl

omit [Fintype Θ] in
/-- **Add-one-in and leave-one-out at once.** The posterior given every unit but data unit `i` is
the `n`-unit posterior reweighted by `g(z_new | θ) / g(z_i | θ)`. -/
theorem postUExcept_some {ι₀ : Type*} [Fintype ι₀] [DecidableEq ι₀] (z : Option ι₀ → Z) (i : ι₀)
    (θ : Θ) (hg : integrated q f θ (z (some i)) ≠ 0) :
    postUExcept prior q f z (some i) θ = postU prior q f (fun j => z (some j)) θ *
      (integrated q f θ (z none) / integrated q f θ (z (some i))) := by
  rw [← mul_div_assoc, eq_div_iff hg, postUExcept_mul, postU_addOneIn]

/-- **The candidate's Hoff score is the plain `n`-unit predictive**, `E_n[g(z_new | θ)]`. -/
theorem hoffScore_none {ι₀ : Type*} [Fintype ι₀] [DecidableEq ι₀] (z : Option ι₀ → Z) :
    hoffScore prior q f z none =
      ∑ θ, post prior q f (fun j => z (some j)) θ * integrated q f θ (z none) := by
  simp only [hoffScore, postUExcept_none, post, div_mul_eq_mul_div, ← sum_div]

/-- `B_i = E_n[g(z_new | θ) / g(z_i | θ)]`, the mean of data unit `i`'s Hoff weight under the
`n`-unit posterior. -/
noncomputable def hoffMeanWeight {ι₀ : Type*} [Fintype ι₀] (z : Option ι₀ → Z) (i : ι₀) : ℝ :=
  ∑ θ, post prior q f (fun j => z (some j)) θ *
    (integrated q f θ (z none) / integrated q f θ (z (some i)))

/-- **Each data unit's Hoff score is the candidate's over its mean weight**: `c_i = c_{n+1} / B_i`.
From `S` draws the estimator is `Â / B̂_i` and `Â`, with `Â` common to every unit, so on the log
scale the bounds slice 10 carries are on `log B̂_i` alone and the candidate's is zero
(`coverage_close_shift`). -/
theorem hoffScore_some {ι₀ : Type*} [Fintype ι₀] [DecidableEq ι₀] (z : Option ι₀ → Z) (i : ι₀)
    (hg : ∀ θ, integrated q f θ (z (some i)) ≠ 0)
    (hN : ∑ θ, postU prior q f (fun j => z (some j)) θ ≠ 0) :
    hoffScore prior q f z (some i) = hoffScore prior q f z none / hoffMeanWeight prior q f z i := by
  rw [hoffScore_eq, hoffScore_none prior q f z, hoffMeanWeight]
  have hnum : ∑ θ, postU prior q f z θ =
      ∑ θ, postU prior q f (fun j => z (some j)) θ * integrated q f θ (z none) :=
    sum_congr rfl fun θ _ => postU_addOneIn prior q f z θ
  have hden : ∑ θ, postUExcept prior q f z (some i) θ =
      ∑ θ, postU prior q f (fun j => z (some j)) θ *
        (integrated q f θ (z none) / integrated q f θ (z (some i))) :=
    sum_congr rfl fun θ _ => postUExcept_some prior q f z i θ (hg θ)
  rw [hnum, hden]
  simp only [post, div_mul_eq_mul_div, ← sum_div]
  rw [div_div_div_cancel_right₀ hN]

/-- **The candidate's accept rule needs only the mean weights**: at level `k` it is accepted iff at
least `k` data units have `B_i ≥ 1`. -/
theorem hoff_accepts_iff {ι₀ : Type*} [Fintype ι₀] [DecidableEq ι₀] (z : Option ι₀ → Z)
    (hg : ∀ θ i, integrated q f θ (z (some i)) ≠ 0)
    (hN : ∑ θ, postU prior q f (fun j => z (some j)) θ ≠ 0)
    (hA : 0 < hoffScore prior q f z none) (hB : ∀ i, 0 < hoffMeanWeight prior q f z i) (k : ℕ) :
    Accepts k univ (hoffScore prior q f z) none ↔
      k ≤ (univ.filter fun i : ι₀ => 1 ≤ hoffMeanWeight prior q f z i).card := by
  have hrank : rankLe univ (hoffScore prior q f z) none =
      1 + (univ.filter fun i : ι₀ => 1 ≤ hoffMeanWeight prior q f z i).card := by
    unfold rankLe
    rw [card_filter, Fintype.sum_option, if_pos le_rfl, ← card_filter]
    congr 1
    refine congrArg card (filter_congr fun i _ => ?_)
    rw [hoffScore_some prior q f z i (fun θ => hg θ i) hN, div_le_iff₀ (hB i),
      le_mul_iff_one_le_right hA]
  unfold Accepts
  rw [hrank]
  constructor <;> intro h <;> omega

/-! #### Hoff's score under latents -/

/-- The unnormalized joint posterior of `θ` and every latent given every unit but `i`: unit `i`'s
latent keeps its prior. -/
noncomputable def jointUExcept (z : ι → Z) (i : ι) (θ : Θ) (v : ι → V) : ℝ :=
  prior θ * ∏ j, (q θ (v j) * if j = i then 1 else f θ (v j) (z j))

/-- **Hoff's score read by the joint posterior given the others**: `E[f(z_i | θ, v_i)]` under the
joint posterior of `θ` and every latent given every unit but `i`. -/
noncomputable def hoffJointScore (z : ι → Z) (i : ι) : ℝ :=
  (∑ θ, ∑ v, jointUExcept prior q f z i θ v * f θ (v i) (z i)) /
    ∑ θ, ∑ v, jointUExcept prior q f z i θ v

/-- **Given the others, unit `i`'s latent is at its prior**, so the joint-posterior reading of
Hoff's score is the integrated one: the latent posterior adds nothing to a leave-one-out score. -/
theorem hoffJointScore_eq (hq1 : ∀ θ, ∑ v, q θ v = 1) (z : ι → Z) (i : ι) :
    hoffJointScore prior q f z i = hoffScore prior q f z i := by
  have ha : ∀ θ, ∀ j ∈ univ.erase i,
      ∑ x, q θ x * (if j = i then 1 else f θ x (z j)) = integrated q f θ (z j) :=
    fun θ j hj => by simp only [if_neg (ne_of_mem_erase hj), integrated]
  have hnum : ∀ θ, ∑ v, jointUExcept prior q f z i θ v * f θ (v i) (z i) =
      postUExcept prior q f z i θ * integrated q f θ (z i) := by
    intro θ
    have h := sum_prod_mul_coord (fun j x => q θ x * (if j = i then 1 else f θ x (z j)))
      (fun x => f θ x (z i)) i
    have hl : ∑ v, jointUExcept prior q f z i θ v * f θ (v i) (z i) =
        prior θ * ∑ v : ι → V,
          (∏ j, q θ (v j) * (if j = i then 1 else f θ (v j) (z j))) * f θ (v i) (z i) := by
      rw [mul_sum]
      exact sum_congr rfl fun v _ => by rw [jointUExcept, mul_assoc]
    rw [hl, h, prod_congr rfl (ha θ)]
    simp only [↓reduceIte, mul_one, postUExcept, integrated]
    ring
  have hden : ∀ θ, ∑ v, jointUExcept prior q f z i θ v = postUExcept prior q f z i θ := by
    intro θ
    have h := sum_prod_mul_coord (fun j x => q θ x * (if j = i then 1 else f θ x (z j)))
      (fun _ => (1 : ℝ)) i
    have hl : ∑ v, jointUExcept prior q f z i θ v = prior θ * ∑ v : ι → V,
        (∏ j, q θ (v j) * (if j = i then 1 else f θ (v j) (z j))) * 1 := by
      rw [mul_sum]
      exact sum_congr rfl fun v _ => by rw [jointUExcept, mul_one]
    rw [hl, h, prod_congr rfl (ha θ)]
    simp only [↓reduceIte, mul_one, hq1, one_mul, postUExcept]
  rw [hoffJointScore, hoffScore, sum_congr rfl fun θ _ => hnum θ,
    sum_congr rfl fun θ _ => hden θ]

/-- **Hoff's score with each unit's latent conditioned on its own value**: `θ` given the other
units, unit `i`'s latent given `θ` and `z_i`, `Σ_θ p(θ | z_{-i}) Σ_v q f(z_i)² / g(z_i)`. This is
what reading the retained latents with Hoff's weights computes (`hoffCondScore_some_draws`). -/
noncomputable def hoffCondScore (z : ι → Z) (i : ι) : ℝ :=
  (∑ θ, postUExcept prior q f z i θ *
      ((∑ v, q θ v * f θ v (z i) ^ 2) / integrated q f θ (z i))) /
    ∑ θ, postUExcept prior q f z i θ

theorem hoffCondScore_comp (z : ι → Z) (σ : Equiv.Perm ι) (i : ι) :
    hoffCondScore prior q f (z ∘ σ) i = hoffCondScore prior q f z (σ i) := by
  simp only [hoffCondScore, postUExcept_comp, Function.comp_apply]

/-- **The conditioned reading stays equivariant**: every unit reads its own latent alike. -/
theorem hoffCondScore_equivariant (C : Finset ι) : Equivariant C (hoffCondScore prior q f) :=
  fun σ _ z i => hoffCondScore_comp prior q f z σ i

omit [Fintype Θ] in
/-- `g ≤ Σ_v q f² / g`: the latent's conditional mean of `f` is at least its prior mean, by the
conditional variance. -/
theorem integrated_le_sq_div (hq : ∀ θ v, 0 ≤ q θ v) (hq1 : ∀ θ, ∑ v, q θ v = 1) (θ : Θ) (y : Z)
    (hg : 0 < integrated q f θ y) :
    integrated q f θ y ≤ (∑ v, q θ v * f θ v y ^ 2) / integrated q f θ y := by
  rw [le_div_iff₀ hg]
  have hvar : 0 ≤ ∑ v, q θ v * (f θ v y - integrated q f θ y) ^ 2 :=
    sum_nonneg fun v _ => mul_nonneg (hq θ v) (sq_nonneg _)
  have hterm : ∀ v, q θ v * (f θ v y - integrated q f θ y) ^ 2 =
      q θ v * f θ v y ^ 2 - 2 * integrated q f θ y * (q θ v * f θ v y) +
        integrated q f θ y ^ 2 * q θ v := fun v => by ring
  have hsum : ∑ v, q θ v * (f θ v y - integrated q f θ y) ^ 2 =
      ∑ v, q θ v * f θ v y ^ 2 - 2 * integrated q f θ y * integrated q f θ y +
        integrated q f θ y ^ 2 * 1 := by
    rw [sum_congr rfl fun v _ => hterm v, sum_add_distrib, sum_sub_distrib, ← mul_sum, ← mul_sum,
      hq1]
    rfl
  rw [hsum] at hvar
  nlinarith [hvar]

omit [Fintype Θ] in
theorem postUExcept_nonneg (hprior : ∀ θ, 0 ≤ prior θ) (hq : ∀ θ v, 0 ≤ q θ v)
    (hf : ∀ θ v y, 0 ≤ f θ v y) (z : ι → Z) (i : ι) (θ : Θ) : 0 ≤ postUExcept prior q f z i θ :=
  mul_nonneg (hprior θ) (prod_nonneg fun j _ =>
    sum_nonneg fun v _ => mul_nonneg (hq θ v) (hf θ v (z j)))

/-- **The conditioned reading sits above the integrated one**, unit by unit. -/
theorem hoffScore_le_hoffCondScore (hprior : ∀ θ, 0 ≤ prior θ) (hq : ∀ θ v, 0 ≤ q θ v)
    (hf : ∀ θ v y, 0 ≤ f θ v y) (hq1 : ∀ θ, ∑ v, q θ v = 1) (z : ι → Z) (i : ι)
    (hg : ∀ θ, 0 < integrated q f θ (z i)) :
    hoffScore prior q f z i ≤ hoffCondScore prior q f z i := by
  refine div_le_div_of_nonneg_right (sum_le_sum fun θ _ => mul_le_mul_of_nonneg_left
    (integrated_le_sq_div q f hq hq1 θ (z i) (hg θ))
    (postUExcept_nonneg prior q f hprior hq hf z i θ)) ?_
  exact sum_nonneg fun θ _ => postUExcept_nonneg prior q f hprior hq hf z i θ

omit [Fintype Θ] in
/-- Every latent but unit `i`'s integrates out of the joint posterior given every unit, leaving
`postU · Σ_v q f(z_i)² / g(z_i)`. -/
theorem sum_jointU_mul_self (z : ι → Z) (i : ι) (θ : Θ)
    (hg : integrated q f θ (z i) ≠ 0) :
    ∑ v, jointU prior q f z θ v * f θ (v i) (z i) = postU prior q f z θ *
      ((∑ v, q θ v * f θ v (z i) ^ 2) / integrated q f θ (z i)) := by
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
  simp only [integrated] at hg ⊢
  field_simp

/-- **The conditioned Hoff score of a data unit from the `n`-unit draws**: the retained latents
read with Hoff's weight `g(z_new | θ) / g(z_i | θ)`. -/
theorem hoffCondScore_some_draws {ι₀ : Type*} [Fintype ι₀] [DecidableEq ι₀] (z : Option ι₀ → Z)
    (i : ι₀) (hg : ∀ θ, integrated q f θ (z (some i)) ≠ 0) :
    hoffCondScore prior q f z (some i) =
      (∑ θ, ∑ v₀ : ι₀ → V, jointU prior q f (fun j => z (some j)) θ v₀ *
          (integrated q f θ (z none) / integrated q f θ (z (some i))) * f θ (v₀ i) (z (some i))) /
        ∑ θ, postU prior q f (fun j => z (some j)) θ *
          (integrated q f θ (z none) / integrated q f θ (z (some i))) := by
  rw [hoffCondScore]
  congr 1
  · refine sum_congr rfl fun θ _ => ?_
    have h := sum_jointU_mul_self prior q f (fun j => z (some j)) i θ (hg θ)
    beta_reduce at h
    rw [postUExcept_some prior q f z i θ (hg θ)]
    calc postU prior q f (fun j => z (some j)) θ *
          (integrated q f θ (z none) / integrated q f θ (z (some i))) *
          ((∑ v, q θ v * f θ v (z (some i)) ^ 2) / integrated q f θ (z (some i)))
        = (integrated q f θ (z none) / integrated q f θ (z (some i))) *
            ∑ v₀ : ι₀ → V, jointU prior q f (fun j => z (some j)) θ v₀ *
              f θ (v₀ i) (z (some i)) := by rw [h]; ring
      _ = _ := by
        rw [mul_sum]
        exact sum_congr rfl fun v₀ _ => by ring
  · exact sum_congr rfl fun θ _ => postUExcept_some prior q f z i θ (hg θ)

omit [Fintype Θ] in
/-- The mean of `f(z_new | θ, v')` under the candidate latent's conditional `q f / g`. -/
theorem sum_cond_mul_self (θ : Θ) (y : Z) :
    ∑ v, q θ v * f θ v y / integrated q f θ y * f θ v y =
      (∑ v, q θ v * f θ v y ^ 2) / integrated q f θ y := by
  rw [sum_div]
  exact sum_congr rfl fun v _ => by ring

/-- **The conditioned Hoff score of the candidate from the `n`-unit draws**: unweighted, with the
candidate's latent drawn from `p(v' | θ, z_new)` for each draw. -/
theorem hoffCondScore_none_draws {ι₀ : Type*} [Fintype ι₀] [DecidableEq ι₀]
    (z : Option ι₀ → Z) :
    hoffCondScore prior q f z none =
      (∑ θ, ∑ v₀ : ι₀ → V, jointU prior q f (fun j => z (some j)) θ v₀ *
          ∑ v', q θ v' * f θ v' (z none) / integrated q f θ (z none) * f θ v' (z none)) /
        ∑ θ, postU prior q f (fun j => z (some j)) θ := by
  rw [hoffCondScore]
  simp only [postUExcept_none, sum_cond_mul_self, ← sum_mul, sum_jointU]

/-- **The mixed Hoff reading**: data units by their conditioned latents, the candidate
integrated. -/
noncomputable def hoffMixedScore {ι₀ : Type*} [Fintype ι₀] [DecidableEq ι₀] (z : Option ι₀ → Z)
    (o : Option ι₀) : ℝ :=
  o.elim (hoffScore prior q f z none) fun j => hoffCondScore prior q f z (some j)

/-- **Reading the units differently breaks Hoff's score as it breaks CB-4's**: the mixed reading
accepts no candidate the conditioned reading refuses. -/
theorem accepts_hoffCond_of_mixed {ι₀ : Type*} [Fintype ι₀] [DecidableEq ι₀]
    (z : Option ι₀ → Z) (hle : hoffScore prior q f z none ≤ hoffCondScore prior q f z none)
    {k : ℕ} (h : Accepts k univ (hoffMixedScore prior q f z) none) :
    Accepts k univ (hoffCondScore prior q f z) none := by
  refine lt_of_lt_of_le h (card_le_card fun j hj => ?_)
  rw [mem_filter] at hj ⊢
  refine ⟨hj.1, ?_⟩
  cases j with
  | none => exact le_rfl
  | some j => exact le_trans hj.2 hle

/-- **The mixed Hoff reading, the counterexample**: CB-4(d)'s model (one `θ`, a fair latent bit,
`y` the bit with probability `3/4`, two units). The data unit scores `5/8` and the candidate
`1/2`, whatever the data, so the candidate is never accepted at `k = 1`: coverage `0` against
CB-1's `1/2`. -/
theorem hoffMixed_counterexample (z : Option Unit → Bool) :
    ¬ Accepts 1 univ (hoffMixedScore (fun _ => (1 : ℝ)) exQ exF z) none := by
  have hg : ∀ y, integrated exQ exF () y = 1 / 2 := by
    intro y
    cases y <;> simp [integrated, exQ, exF] <;> norm_num
  have hE : ∀ i, postUExcept (fun _ => (1 : ℝ)) exQ exF z i () ≠ 0 := by
    intro i
    rw [postUExcept, one_mul]
    exact prod_ne_zero_iff.2 fun j _ => by rw [hg]; norm_num
  have hU : ∀ F : Unit → ℝ, ∑ θ, F θ = F () := fun F => Fintype.sum_unique F
  have hint : hoffMixedScore (fun _ => (1 : ℝ)) exQ exF z none = 1 / 2 := by
    change hoffScore (fun _ => (1 : ℝ)) exQ exF z none = 1 / 2
    rw [hoffScore, hU, hU, mul_div_cancel_left₀ _ (hE none), hg]
  have hlat : ∀ u, hoffMixedScore (fun _ => (1 : ℝ)) exQ exF z (some u) = 5 / 8 := by
    intro u
    change hoffCondScore (fun _ => (1 : ℝ)) exQ exF z (some u) = 5 / 8
    rw [hoffCondScore, hU, hU, mul_div_cancel_left₀ _ (hE (some u)), hg]
    simp only [exQ, exF, Fintype.sum_bool]
    cases z (some u) <;> norm_num
  have hfilt : (univ.filter fun j => hoffMixedScore (fun _ => (1 : ℝ)) exQ exF z j ≤
      hoffMixedScore (fun _ => (1 : ℝ)) exQ exF z none) = {none} := by
    ext j
    cases j with
    | none => simp
    | some u =>
      simp only [mem_filter, mem_univ, true_and, mem_singleton, reduceCtorEq, iff_false, not_le,
        hint, hlat]
      norm_num
  simp only [Accepts, rankLe, hfilt, card_singleton, lt_irrefl, not_false_eq_true]

/-! #### CB-4(c) read from draws, the candidate's latent drawn per draw -/

/-- **CB-4(c) of a data unit from the `n`-unit draws**: the retained latents, read with the
add-one-in weight `g(z_new | θ)`. -/
theorem latentScore_some_draws {ι₀ : Type*} [Fintype ι₀] [DecidableEq ι₀] (z : Option ι₀ → Z)
    (i : ι₀) :
    latentScore prior q f z (some i) =
      (∑ θ, ∑ v₀ : ι₀ → V, jointU prior q f (fun j => z (some j)) θ v₀ *
          integrated q f θ (z none) * f θ (v₀ i) (z (some i))) /
        ∑ θ, postU prior q f (fun j => z (some j)) θ * integrated q f θ (z none) := by
  unfold latentScore
  congr 1
  · refine sum_congr rfl fun θ _ => ?_
    rw [sum_jointU_addOneIn prior q f z θ (fun v => f θ (v (some i)) (z (some i)))]
    refine sum_congr rfl fun v₀ _ => ?_
    simp only [Option.elim_some]
    rw [← sum_mul, integrated]
    ring
  · exact sum_congr rfl fun θ _ => postU_addOneIn prior q f z θ

/-- **CB-4(c) of the candidate from the `n`-unit draws**: the add-one-in weight, and the
candidate's latent drawn from `p(v' | θ, z_new)` for each draw. So the estimator's expectation
over the drawn latents (`drawnLatent_mean_var`) gives (c)'s score for every unit, the candidate
included. -/
theorem latentScore_none_draws {ι₀ : Type*} [Fintype ι₀] [DecidableEq ι₀] (z : Option ι₀ → Z)
    (hg : ∀ θ, integrated q f θ (z none) ≠ 0) :
    latentScore prior q f z none =
      (∑ θ, ∑ v₀ : ι₀ → V, jointU prior q f (fun j => z (some j)) θ v₀ *
          integrated q f θ (z none) *
            ∑ v', q θ v' * f θ v' (z none) / integrated q f θ (z none) * f θ v' (z none)) /
        ∑ θ, postU prior q f (fun j => z (some j)) θ * integrated q f θ (z none) := by
  unfold latentScore
  congr 1
  · refine sum_congr rfl fun θ _ => ?_
    rw [sum_jointU_addOneIn prior q f z θ (fun v => f θ (v none) (z none))]
    refine sum_congr rfl fun v₀ _ => ?_
    simp only [Option.elim_none]
    rw [sum_cond_mul_self, mul_assoc, ← mul_div_assoc, mul_div_cancel_left₀ _ (hg θ)]
    congr 1
    exact sum_congr rfl fun v _ => by ring
  · exact sum_congr rfl fun θ _ => postU_addOneIn prior q f z θ

end Hoff

section DrawnLatents

variable {S V : Type*} [Fintype S] [DecidableEq S] [Fintype V]

/-- Independent draws factor: the sum over every draw's value of a product of per-draw factors. -/
theorem sum_pi_prod_mul_prod (a b : S → V → ℝ) :
    ∑ v : S → V, (∏ j, a j (v j)) * ∏ j, b j (v j) = ∏ j, ∑ x, a j x * b j x := by
  have h := Finset.prod_univ_sum (t := fun _ : S => (univ : Finset V))
    (f := fun j x => a j x * b j x)
  rw [Fintype.piFinset_univ] at h
  rw [h]
  exact sum_congr rfl fun v _ => prod_mul_distrib.symm

omit [Fintype S] [Fintype V] in
/-- The factor reading `H` at draw `s`, `H'` at draw `s'`, and nothing elsewhere. -/
def twoPoint (s s' : S) (H H' : V → ℝ) (j : S) (x : V) : ℝ :=
  if j = s then H x else if j = s' then H' x else 1

omit [Fintype V] in
theorem prod_twoPoint {s s' : S} (hss : s ≠ s') (H H' : V → ℝ) (v : S → V) :
    ∏ j, twoPoint s s' H H' j (v j) = H (v s) * H' (v s') := by
  rw [Fintype.prod_eq_mul s s' hss fun j hj => by simp only [twoPoint, if_neg hj.1, if_neg hj.2]]
  simp only [twoPoint, ↓reduceIte, if_neg hss.symm]

/-- Two distinct draws are independent. -/
theorem sum_pi_prod_mul_two (p : S → V → ℝ) (hp1 : ∀ s, ∑ v, p s v = 1) {s s' : S}
    (hss : s ≠ s') (H H' : V → ℝ) :
    ∑ v : S → V, (∏ j, p j (v j)) * (H (v s) * H' (v s')) =
      (∑ x, p s x * H x) * ∑ x, p s' x * H' x := by
  have e : ∀ v : S → V, (∏ j, p j (v j)) * (H (v s) * H' (v s')) =
      (∏ j, p j (v j)) * ∏ j, twoPoint s s' H H' j (v j) := fun v => by rw [prod_twoPoint hss]
  rw [sum_congr rfl fun v _ => e v, sum_pi_prod_mul_prod,
    Fintype.prod_eq_mul s s' hss fun j hj => by
      simp only [twoPoint, if_neg hj.1, if_neg hj.2, mul_one, hp1]]
  simp only [twoPoint, ↓reduceIte, if_neg hss.symm]

/-- The mean of a sum of per-draw readings, each draw's latent drawn independently. -/
theorem sum_pi_prod_mul_sum (p : S → V → ℝ) (hp1 : ∀ s, ∑ v, p s v = 1) (h : S → V → ℝ) :
    ∑ v : S → V, (∏ j, p j (v j)) * ∑ s, h s (v s) = ∑ s, ∑ x, p s x * h s x := by
  simp only [mul_sum]
  refine sum_comm.trans (sum_congr rfl fun s _ => ?_)
  rw [sum_prod_mul_coord p (h s) s, prod_congr rfl fun j _ => hp1 j, prod_const_one, mul_one]

/-- The variance of a sum of centred per-draw readings is the sum of their variances. -/
theorem sum_pi_prod_mul_sq (p : S → V → ℝ) (hp1 : ∀ s, ∑ v, p s v = 1) (d : S → V → ℝ)
    (hd : ∀ s, ∑ x, p s x * d s x = 0) :
    ∑ v : S → V, (∏ j, p j (v j)) * (∑ s, d s (v s)) ^ 2 = ∑ s, ∑ x, p s x * d s x ^ 2 := by
  have hcross : ∀ s s', ∑ v : S → V, (∏ j, p j (v j)) * (d s (v s) * d s' (v s')) =
      if s = s' then ∑ x, p s x * d s x ^ 2 else 0 := by
    intro s s'
    split_ifs with h
    · rw [← h]
      have hc := sum_prod_mul_coord p (fun x => d s x ^ 2) s
      beta_reduce at hc
      rw [prod_congr rfl fun j _ => hp1 j, prod_const_one, mul_one] at hc
      rw [← hc]
      exact sum_congr rfl fun v _ => by ring
    · rw [sum_pi_prod_mul_two p hp1 h, hd, zero_mul]
  calc ∑ v : S → V, (∏ j, p j (v j)) * (∑ s, d s (v s)) ^ 2
      = ∑ v : S → V, ∑ s, ∑ s', (∏ j, p j (v j)) * (d s (v s) * d s' (v s')) := by
        refine sum_congr rfl fun v _ => ?_
        rw [sq, sum_mul_sum, mul_sum]
        exact sum_congr rfl fun s _ => by rw [mul_sum]
    _ = ∑ s, ∑ s', ∑ v : S → V, (∏ j, p j (v j)) * (d s (v s) * d s' (v s')) := by
        exact sum_comm.trans (sum_congr rfl fun s _ => sum_comm)
    _ = ∑ s, ∑ x, p s x * d s x ^ 2 := by
        refine sum_congr rfl fun s _ => ?_
        rw [sum_congr rfl fun s' _ => hcross s s', sum_ite_eq, if_pos (mem_univ s)]

/-- **What the drawn latent adds.** Draws `s` carry normalized weights `a s` (the add-one-in
weights over their sum), the candidate's latent at draw `s` is drawn from `p s` (its conditional
`p(v' | θ_s, z_new)`), and `h s v` is its reading `f(z_new | θ_s, v)`. Over the drawn latents the
estimate `Σ_s a_s h_s(v_s)` has mean `Σ_s a_s E_s h_s`, the Rao–Blackwell reading, and variance
`Σ_s a_s² Var_s h_s`: the term the candidate's bound must carry on top of the draws' own error,
and that reading the conditional mean in closed form removes. -/
theorem drawnLatent_mean_var (p : S → V → ℝ) (hp1 : ∀ s, ∑ v, p s v = 1) (a : S → ℝ)
    (h : S → V → ℝ) :
    ∑ v : S → V, (∏ j, p j (v j)) * ∑ s, a s * h s (v s) = ∑ s, a s * ∑ x, p s x * h s x ∧
    ∑ v : S → V, (∏ j, p j (v j)) * (∑ s, a s * (h s (v s) - ∑ x, p s x * h s x)) ^ 2 =
      ∑ s, a s ^ 2 * ∑ x, p s x * (h s x - ∑ x', p s x' * h s x') ^ 2 := by
  constructor
  · rw [sum_pi_prod_mul_sum p hp1 (fun s x => a s * h s x)]
    exact sum_congr rfl fun s _ => by rw [mul_sum]; exact sum_congr rfl fun x _ => by ring
  · have hd : ∀ s, ∑ x, p s x * (a s * (h s x - ∑ x', p s x' * h s x')) = 0 := by
      intro s
      have e : ∀ x, p s x * (a s * (h s x - ∑ x', p s x' * h s x')) =
          a s * (p s x * h s x) - (a s * ∑ x', p s x' * h s x') * p s x := fun x => by ring
      rw [sum_congr rfl fun x _ => e x, sum_sub_distrib, ← mul_sum, ← mul_sum, hp1]
      ring
    rw [sum_pi_prod_mul_sq p hp1 (fun s x => a s * (h s x - ∑ x', p s x' * h s x')) hd]
    exact sum_congr rfl fun s _ => by rw [mul_sum]; exact sum_congr rfl fun x _ => by ring

end DrawnLatents

section HoffOptimal

variable {ι Z : Type*} [Fintype ι] [DecidableEq ι] [Fintype Z]

/-- The others' mass: the law of every unit but `i`, `P_{-i}(z) = Σ_y P(z with i ↦ y)`. -/
noncomputable def othersMass (P : (ι → Z) → ℝ) (i : ι) (z : ι → Z) : ℝ :=
  ∑ y, P (Function.update z i y)

/-- **Hoff's score for any law `P` of the units**: the `μ`-density of unit `i`'s value given the
others, `c_i = P(z) / (μ(z_i) P_{-i}(z))`. -/
noncomputable def hoffP (P : (ι → Z) → ℝ) (μ : Z → ℝ) (z : ι → Z) (i : ι) : ℝ :=
  P z / (μ (z i) * othersMass P i z)

/-- **A procedure's expected size** under the law `P`. `φ z` is the chance that the candidate's
value `z t` is accepted given the others, and the size counts accepted values by `μ`. Summing over
every `z` counts each value of the others once per candidate value, as `P_{-t}` does not read
`z t`. -/
noncomputable def expectedSize (P : (ι → Z) → ℝ) (μ : Z → ℝ) (t : ι) (φ : (ι → Z) → ℝ) : ℝ :=
  ∑ z, othersMass P t z * μ (z t) * φ z

omit [Fintype ι] [Fintype Z] in
theorem update_comp_perm (z : ι → Z) (σ : Equiv.Perm ι) (i : ι) (y : Z) :
    Function.update (z ∘ σ) i y = Function.update z (σ i) y ∘ σ := by
  funext j
  simp only [Function.comp_apply, Function.update_apply]
  by_cases h : j = i
  · simp [h]
  · have : σ j ≠ σ i := fun h' => h (σ.injective h')
    simp [h, this]

omit [Fintype ι] in
theorem othersMass_comp {P : (ι → Z) → ℝ} (hP : ∀ σ : Equiv.Perm ι, ∀ z, P (z ∘ σ) = P z)
    (z : ι → Z) (σ : Equiv.Perm ι) (i : ι) :
    othersMass P i (z ∘ σ) = othersMass P (σ i) z := by
  unfold othersMass
  exact sum_congr rfl fun y _ => by rw [update_comp_perm, hP]

omit [Fintype ι] in
/-- Hoff's score is equivariant under any exchangeable law. -/
theorem hoffP_comp {P : (ι → Z) → ℝ} (hP : ∀ σ : Equiv.Perm ι, ∀ z, P (z ∘ σ) = P z)
    (μ : Z → ℝ) (z : ι → Z) (σ : Equiv.Perm ι) (i : ι) :
    hoffP P μ (z ∘ σ) i = hoffP P μ z (σ i) := by
  simp only [hoffP, hP, Function.comp_apply, othersMass_comp hP]

omit [Fintype ι] in
theorem othersMass_pos {P : (ι → Z) → ℝ} (hpos : ∀ z, 0 < P z) (i : ι) (z : ι → Z) :
    0 < othersMass P i z :=
  sum_pos (fun _ _ => hpos _) ⟨z i, mem_univ _⟩

omit [Fintype ι] in
theorem hoffP_pos {P : (ι → Z) → ℝ} (hpos : ∀ z, 0 < P z) {μ : Z → ℝ} (hμ : ∀ y, 0 < μ y)
    (z : ι → Z) (i : ι) : 0 < hoffP P μ z i :=
  div_pos (hpos z) (mul_pos (hμ _) (othersMass_pos hpos i z))

theorem expectedSize_eq {P : (ι → Z) → ℝ} (hpos : ∀ z, 0 < P z) (μ : Z → ℝ) (t : ι)
    (χ : (ι → Z) → ℝ) :
    expectedSize P μ t χ = ∑ z, P z * (χ z / hoffP P μ z t) := by
  refine sum_congr rfl fun z _ => ?_
  rw [hoffP, div_div_eq_mul_div, ← mul_div_assoc, mul_div_cancel_left₀ _ (hpos z).ne']
  ring

theorem sum_comp_perm (F : (ι → Z) → ℝ) (σ : Equiv.Perm ι) :
    ∑ z, F (z ∘ σ) = ∑ z, F z := by
  let e : (ι → Z) ≃ (ι → Z) :=
    ⟨fun z => z ∘ σ, fun z => z ∘ σ.symm, fun z => by funext i; simp, fun z => by funext i; simp⟩
  exact Fintype.sum_equiv e (fun z => F (z ∘ σ)) F fun z => rfl

/-- **The expected size, averaged over the orderings.** By exchangeability each data set's share
of the size is its orderings' readings weighted by `1 / c` at the unit the candidate lands on. -/
theorem sum_perm_expectedSize {P : (ι → Z) → ℝ} (hP : ∀ σ : Equiv.Perm ι, ∀ z, P (z ∘ σ) = P z)
    (hpos : ∀ z, 0 < P z) (μ : Z → ℝ) (t : ι) (χ : (ι → Z) → ℝ) :
    ∑ _σ : Equiv.Perm ι, expectedSize P μ t χ =
      ∑ z, P z * ∑ σ : Equiv.Perm ι, χ (z ∘ σ) / hoffP P μ z (σ t) := by
  rw [expectedSize_eq hpos μ t χ]
  calc ∑ _σ : Equiv.Perm ι, ∑ z, P z * (χ z / hoffP P μ z t)
      = ∑ σ : Equiv.Perm ι, ∑ z, P (z ∘ σ) * (χ (z ∘ σ) / hoffP P μ (z ∘ σ) t) :=
        sum_congr rfl fun σ _ => (sum_comp_perm (fun z => P z * (χ z / hoffP P μ z t)) σ).symm
    _ = ∑ z, ∑ σ : Equiv.Perm ι, P z * (χ (z ∘ σ) / hoffP P μ z (σ t)) := by
        rw [sum_comm]
        exact sum_congr rfl fun z _ => sum_congr rfl fun σ _ => by rw [hP, hoffP_comp hP]
    _ = _ := sum_congr rfl fun z _ => by rw [mul_sum]

omit [Fintype ι] [DecidableEq ι] [Fintype Z] in
/-- **The bag's linear program.** Readings `a` and `b` in `[0, 1]`, with `b` a threshold in the
weights `w ≥ 0` (wherever `b x > 0` and `b y < 1`, `w x ≤ w y`) and no more total than `a`: then
`b`'s weighted total is no more than `a`'s. -/
theorem sum_mul_le_of_threshold {α : Type*} (s : Finset α) (a b w : α → ℝ)
    (ha0 : ∀ x ∈ s, 0 ≤ a x) (ha1 : ∀ x ∈ s, a x ≤ 1) (hb0 : ∀ x ∈ s, 0 ≤ b x)
    (hb1 : ∀ x ∈ s, b x ≤ 1) (hw0 : ∀ x ∈ s, 0 ≤ w x)
    (hsep : ∀ x ∈ s, ∀ y ∈ s, 0 < b x → b y < 1 → w x ≤ w y)
    (hcount : ∑ x ∈ s, b x ≤ ∑ x ∈ s, a x) :
    ∑ x ∈ s, b x * w x ≤ ∑ x ∈ s, a x * w x := by
  classical
  rcases (s.filter fun x => 0 < b x).eq_empty_or_nonempty with hT | hT
  · have hzero : ∀ x ∈ s, b x = 0 := fun x hx => by
      by_contra hne
      have hmem : x ∈ s.filter fun x => 0 < b x :=
        mem_filter.2 ⟨hx, lt_of_le_of_ne (hb0 x hx) (Ne.symm hne)⟩
      rw [hT] at hmem
      simp at hmem
    rw [sum_eq_zero fun x hx => by rw [hzero x hx, zero_mul]]
    exact sum_nonneg fun x hx => mul_nonneg (ha0 x hx) (hw0 x hx)
  · set τ := (s.filter fun x => 0 < b x).sup' hT w
    have hle : ∀ x ∈ s, 0 < b x → w x ≤ τ := fun x hx hb => le_sup' w (mem_filter.2 ⟨hx, hb⟩)
    have hge : ∀ y ∈ s, b y < 1 → τ ≤ w y := fun y hy hb =>
      sup'_le hT w fun x hx => by
        rw [mem_filter] at hx
        exact hsep x hx.1 y hy hx.2 hb
    have hτ : 0 ≤ τ := by
      obtain ⟨x, hx⟩ := id hT
      rw [mem_filter] at hx
      exact (hw0 x hx.1).trans (hle x hx.1 hx.2)
    have key : ∀ x ∈ s, (a x - b x) * τ ≤ (a x - b x) * w x := by
      intro x hx
      rcases le_or_gt (b x) (a x) with hab | hab
      · rcases lt_or_eq_of_le (hb1 x hx) with hb | hb
        · exact mul_le_mul_of_nonneg_left (hge x hx hb) (by linarith)
        · have : a x = b x := le_antisymm (by rw [hb]; exact ha1 x hx) hab
          rw [this, sub_self, zero_mul, zero_mul]
      · exact mul_le_mul_of_nonpos_left (hle x hx (lt_of_le_of_lt (ha0 x hx) hab)) (by linarith)
    have hsum := sum_le_sum key
    rw [← sum_mul] at hsum
    have hnn : 0 ≤ (∑ x ∈ s, (a x - b x)) * τ :=
      mul_nonneg (by rw [sum_sub_distrib]; linarith) hτ
    have h2 : ∑ x ∈ s, (a x - b x) * w x = ∑ x ∈ s, a x * w x - ∑ x ∈ s, b x * w x := by
      rw [← sum_sub_distrib]
      exact sum_congr rfl fun x _ => by ring
    linarith

/-- **Hoff's optimality, finite** (Hoff 2023, Theorem 4.1). Under an exchangeable, positive law
`P` of the units, take a procedure `ψ` that on every data set is a threshold in Hoff's score:
wherever it gives some chance to the candidate landing on one unit and less than certainty on
another, the first unit's Hoff score is at least the second's. CB-1's rule on Hoff's score is one
(`hoff_cb1_expectedSize_le`), and so is its form with ties broken at random. Every procedure `φ`
whose coverage over the orderings of every data set is at least `ψ`'s has an expected size at
least `ψ`'s. `φ` may be any procedure, randomized or not, symmetric in the data or not. -/
theorem hoff_expectedSize_le {P : (ι → Z) → ℝ} (hP : ∀ σ : Equiv.Perm ι, ∀ z, P (z ∘ σ) = P z)
    (hpos : ∀ z, 0 < P z) {μ : Z → ℝ} (hμ : ∀ y, 0 < μ y) (t : ι) (φ ψ : (ι → Z) → ℝ)
    (hφ0 : ∀ z, 0 ≤ φ z) (hφ1 : ∀ z, φ z ≤ 1) (hψ0 : ∀ z, 0 ≤ ψ z) (hψ1 : ∀ z, ψ z ≤ 1)
    (hψ : ∀ z, ∀ σ σ' : Equiv.Perm ι, 0 < ψ (z ∘ σ) → ψ (z ∘ σ') < 1 →
      hoffP P μ z (σ' t) ≤ hoffP P μ z (σ t))
    (hcov : ∀ z, ∑ σ : Equiv.Perm ι, ψ (z ∘ σ) ≤ ∑ σ : Equiv.Perm ι, φ (z ∘ σ)) :
    expectedSize P μ t ψ ≤ expectedSize P μ t φ := by
  have hbag : ∀ z, ∑ σ : Equiv.Perm ι, ψ (z ∘ σ) / hoffP P μ z (σ t) ≤
      ∑ σ : Equiv.Perm ι, φ (z ∘ σ) / hoffP P μ z (σ t) := by
    intro z
    simp only [div_eq_mul_inv]
    refine sum_mul_le_of_threshold (univ : Finset (Equiv.Perm ι)) (fun σ => φ (z ∘ σ))
      (fun σ => ψ (z ∘ σ)) (fun σ => (hoffP P μ z (σ t))⁻¹) (fun σ _ => hφ0 _) (fun σ _ => hφ1 _)
      (fun σ _ => hψ0 _) (fun σ _ => hψ1 _)
      (fun σ _ => (inv_pos.2 (hoffP_pos hpos hμ z (σ t))).le)
      (fun σ _ σ' _ h0 h1 => inv_anti₀ (hoffP_pos hpos hμ z (σ' t)) (hψ z σ σ' h0 h1))
      (hcov z)
  have h := sum_le_sum fun z (_ : z ∈ univ) =>
    mul_le_mul_of_nonneg_left (hbag z) (hpos z).le
  rw [← sum_perm_expectedSize hP hpos μ t ψ, ← sum_perm_expectedSize hP hpos μ t φ,
    sum_const, sum_const, nsmul_eq_mul, nsmul_eq_mul] at h
  exact le_of_mul_le_mul_left h (by exact_mod_cast card_pos.2 ⟨1, mem_univ _⟩)

/-- **CB-1's rule on Hoff's score is smallest** among procedures that cover at least as much on
every data set. With Hoff's scores distinct on a data set that is exactly `1 − k/m`, so then the
comparison is with every procedure valid at `1 − k/m`. -/
theorem hoff_cb1_expectedSize_le {P : (ι → Z) → ℝ}
    (hP : ∀ σ : Equiv.Perm ι, ∀ z, P (z ∘ σ) = P z) (hpos : ∀ z, 0 < P z) {μ : Z → ℝ}
    (hμ : ∀ y, 0 < μ y) (t : ι) (k : ℕ) (φ : (ι → Z) → ℝ) (hφ0 : ∀ z, 0 ≤ φ z)
    (hφ1 : ∀ z, φ z ≤ 1)
    (hcov : ∀ z, (((univ : Finset (Equiv.Perm ι)).filter fun σ : Equiv.Perm ι =>
      Accepts k univ (hoffP P μ (z ∘ σ)) t).card : ℝ) ≤ ∑ σ : Equiv.Perm ι, φ (z ∘ σ)) :
    expectedSize P μ t (fun z => if Accepts k univ (hoffP P μ z) t then 1 else 0) ≤
      expectedSize P μ t φ := by
  have hacc : ∀ z (σ : Equiv.Perm ι),
      Accepts k univ (hoffP P μ (z ∘ σ)) t ↔ Accepts k univ (hoffP P μ z) (σ t) := by
    intro z σ
    have hfun : hoffP P μ (z ∘ σ) = hoffP P μ z ∘ σ := funext fun i => hoffP_comp hP μ z σ i
    rw [hfun]
    unfold Accepts
    rw [rankLe_comp (mem_permsOn_univ σ)]
  refine hoff_expectedSize_le hP hpos hμ t φ _ hφ0 hφ1
    (fun z => by split_ifs <;> norm_num)
    (fun z => by split_ifs <;> norm_num) (fun z σ σ' h0 h1 => ?_) (fun z => ?_)
  · beta_reduce at h0 h1
    have ha : Accepts k univ (hoffP P μ z) (σ t) := by
      rw [← hacc]
      by_contra hn
      rw [if_neg hn] at h0
      exact lt_irrefl _ h0
    have hr : ¬ Accepts k univ (hoffP P μ z) (σ' t) := by
      rw [← hacc]
      intro hn
      rw [if_pos hn] at h1
      exact lt_irrefl _ h1
    by_contra hlt
    rw [not_le] at hlt
    apply hr
    refine lt_of_lt_of_le ha (card_le_card fun j hj => ?_)
    rw [mem_filter] at hj ⊢
    exact ⟨hj.1, hj.2.trans hlt.le⟩
  · rw [sum_boole]
    exact hcov z

end HoffOptimal

section HoffLink

variable {Θ V Z : Type*} [Fintype Θ] [Fintype V] [Fintype Z]
variable (prior : Θ → ℝ) (q : Θ → V → ℝ) (f : Θ → V → Z → ℝ)
variable {ι : Type*} [Fintype ι] [DecidableEq ι]

omit [Fintype Z] in
theorem prod_update_erase (F : Z → ℝ) (z : ι → Z) (i : ι) (y : Z) :
    ∏ j, F (Function.update z i y j) = F y * ∏ j ∈ univ.erase i, F (z j) := by
  rw [← mul_prod_erase univ _ (mem_univ i), Function.update_self]
  congr 1
  exact prod_congr rfl fun j hj => by rw [Function.update_of_ne (ne_of_mem_erase hj)]

/-- **The model's Hoff score is `hoffP` of its law.** With the units' law
`P(z) = Σ_θ π(θ) ∏_j μ(z_j) g(z_j | θ)` and `g` a `μ`-density, `hoffP P μ` is `hoffScore`, so
`hoff_expectedSize_le` is Hoff's optimality for the model's own score. -/
theorem hoffScore_eq_hoffP (μ : Z → ℝ) (hμ : ∀ y, μ y ≠ 0)
    (hg1 : ∀ θ, ∑ y, μ y * integrated q f θ y = 1) (z : ι → Z) (i : ι) :
    hoffScore prior q f z i = hoffP (fun w => ∑ θ, prior θ * ∏ j, (μ (w j) *
      integrated q f θ (w j))) μ z i := by
  have hP : (∑ θ, prior θ * ∏ j, (μ (z j) * integrated q f θ (z j))) =
      (μ (z i) * ∏ j ∈ univ.erase i, μ (z j)) * ∑ θ, postU prior q f z θ := by
    rw [mul_prod_erase univ (fun j => μ (z j)) (mem_univ i), mul_sum]
    refine sum_congr rfl fun θ _ => ?_
    rw [postU, prod_mul_distrib]
    ring
  have hM : othersMass (fun w => ∑ θ, prior θ * ∏ j, (μ (w j) * integrated q f θ (w j))) i z =
      (∏ j ∈ univ.erase i, μ (z j)) * ∑ θ, postUExcept prior q f z i θ := by
    unfold othersMass
    beta_reduce
    rw [sum_comm, mul_sum]
    refine sum_congr rfl fun θ _ => ?_
    calc ∑ y, prior θ * ∏ j, (μ (Function.update z i y j) *
          integrated q f θ (Function.update z i y j))
        = ∑ y, prior θ * ((μ y * integrated q f θ y) *
            ∏ j ∈ univ.erase i, (μ (z j) * integrated q f θ (z j))) :=
          sum_congr rfl fun y _ => by
            rw [prod_update_erase (fun w => μ w * integrated q f θ w) z i y]
      _ = (∏ j ∈ univ.erase i, μ (z j)) * postUExcept prior q f z i θ := by
          rw [← mul_sum, ← sum_mul, hg1, one_mul, postUExcept, prod_mul_distrib]
          ring
  have hne : ∏ j ∈ univ.erase i, μ (z j) ≠ 0 := prod_ne_zero_iff.2 fun j _ => hμ _
  unfold hoffP
  beta_reduce
  rw [hoffScore_eq, hP, hM, ← mul_assoc, mul_div_mul_left _ _ (mul_ne_zero (hμ _) hne)]

end HoffLink

/-! ### Slice 10, the follow-up: reading from the realized ordering

The engine holds one ordering's estimates, the realized one. When that ordering's own test unit
is in `D`, its estimates carry no guarantee, for any unit: every unit's estimate shares the test
unit's add-one-in weights. Such an ordering reads nothing. -/

section ReadingExcluded

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- **A count read from the realized ordering, from above.** `c w ≤ B` is a count of the data set
under the draw `w`, and `ĉ w σ ≥ 0` is what the ordering `σ` reads. The reading is at least `c w`
whenever `σ`'s test unit is kept (outside `D w`) and `σ`'s bounds hold, and a left-out ordering
reads `0`. Then `c` is at most the orderings' average reading, plus `B` times `D`'s share and the
failure share `ν`. -/
theorem le_read_excluding {W : Type*} [Fintype W] {C : Finset ι} {t : ι} (ht : t ∈ C)
    (p : W → ℝ) (hp0 : ∀ w, 0 ≤ p w) (D : W → Finset ι) (hold : W → Equiv.Perm ι → Prop)
    [∀ w, DecidablePred (hold w)] (c : W → ℝ) (ĉ : W → Equiv.Perm ι → ℝ) {B : ℝ} (hB : 0 ≤ B)
    (hcB : ∀ w, c w ≤ B) (hĉ0 : ∀ w σ, 0 ≤ ĉ w σ)
    (hread : ∀ w, ∀ σ ∈ permsOn C, σ t ∉ D w → hold w σ → c w ≤ ĉ w σ) {ν : ℝ}
    (hfail : ∑ w, p w * (((permsOn C).filter fun σ : Equiv.Perm ι =>
        σ t ∉ D w ∧ ¬ hold w σ).card : ℝ) ≤ ν * (permsOn C).card) :
    ∑ w, p w * c w ≤
      ∑ w, p w * ((∑ σ ∈ permsOn C, if σ t ∈ D w then 0 else ĉ w σ) / (permsOn C).card) +
        B * (∑ w, p w * (((C.filter fun u => u ∈ D w).card : ℝ) / C.card) + ν) := by
  have hg : (0 : ℝ) < (permsOn C).card := by exact_mod_cast card_permsOn_pos C
  have hw : ∀ w, c w ≤ (∑ σ ∈ permsOn C, if σ t ∈ D w then 0 else ĉ w σ) / (permsOn C).card +
      B * (((C.filter fun u => u ∈ D w).card : ℝ) / C.card) +
      B * ((((permsOn C).filter fun σ : Equiv.Perm ι => σ t ∉ D w ∧ ¬ hold w σ).card : ℝ) /
        (permsOn C).card) := by
    intro w
    have hσ : ∀ σ ∈ permsOn C, c w ≤ (if σ t ∈ D w then 0 else ĉ w σ) +
        B * (if σ t ∈ D w then 1 else 0) + B * (if σ t ∉ D w ∧ ¬ hold w σ then 1 else 0) := by
      intro σ hσ
      split_ifs with h1 h2 h2
      · linarith [hcB w]
      · linarith [hcB w]
      · linarith [hĉ0 w σ, hcB w]
      · have hh : hold w σ := by
          by_contra hh
          exact h2 ⟨h1, hh⟩
        linarith [hread w σ hσ h1 hh]
    have hsum := sum_le_sum hσ
    rw [sum_const, nsmul_eq_mul, sum_add_distrib, sum_add_distrib, ← mul_sum, ← mul_sum,
      sum_boole, sum_boole] at hsum
    rw [← ratio_permsOn ht (fun u => u ∈ D w)]
    have e : (∑ σ ∈ permsOn C, if σ t ∈ D w then 0 else ĉ w σ) / (permsOn C).card +
        B * ((((permsOn C).filter fun σ : Equiv.Perm ι => σ t ∈ D w).card : ℝ) /
          (permsOn C).card) +
        B * ((((permsOn C).filter fun σ : Equiv.Perm ι => σ t ∉ D w ∧ ¬ hold w σ).card : ℝ) /
          (permsOn C).card) =
        ((∑ σ ∈ permsOn C, if σ t ∈ D w then 0 else ĉ w σ) +
          B * (((permsOn C).filter fun σ : Equiv.Perm ι => σ t ∈ D w).card : ℝ) +
          B * (((permsOn C).filter fun σ : Equiv.Perm ι => σ t ∉ D w ∧ ¬ hold w σ).card : ℝ)) /
            (permsOn C).card := by ring
    rw [e, le_div_iff₀ hg, mul_comm (c w)]
    linarith
  have hsum := sum_le_sum fun w (_ : w ∈ univ) => mul_le_mul_of_nonneg_left (hw w) (hp0 w)
  simp only [mul_add, sum_add_distrib] at hsum
  have hν : ∑ w, p w * ((((permsOn C).filter fun σ : Equiv.Perm ι =>
      σ t ∉ D w ∧ ¬ hold w σ).card : ℝ) / (permsOn C).card) ≤ ν :=
    weighted_share_le C p (fun w σ => σ t ∉ D w ∧ ¬ hold w σ) hfail
  have e1 : ∑ w, p w * (B * (((C.filter fun u => u ∈ D w).card : ℝ) / C.card)) =
      B * ∑ w, p w * (((C.filter fun u => u ∈ D w).card : ℝ) / C.card) := by
    rw [mul_sum]
    exact sum_congr rfl fun w _ => by ring
  have e2 : ∑ w, p w * (B * ((((permsOn C).filter fun σ : Equiv.Perm ι =>
      σ t ∉ D w ∧ ¬ hold w σ).card : ℝ) / (permsOn C).card)) =
      B * ∑ w, p w * ((((permsOn C).filter fun σ : Equiv.Perm ι =>
        σ t ∉ D w ∧ ¬ hold w σ).card : ℝ) / (permsOn C).card) := by
    rw [mul_sum]
    exact sum_congr rfl fun w _ => by ring
  have e3 := mul_le_mul_of_nonneg_left hν hB
  rw [mul_add]
  linarith

/-- **A count read from the realized ordering, from below.** A reading at most `B`, at most the
count `c ≥ 0` whenever the ordering is kept and its bounds hold, and `0` on a left-out ordering:
the orderings' average reading is at most `c` plus `B ν`. No `D` share is charged, since a
left-out ordering reads `0`. -/
theorem read_le_excluding {W : Type*} [Fintype W] {C : Finset ι} {t : ι}
    (p : W → ℝ) (hp0 : ∀ w, 0 ≤ p w) (D : W → Finset ι) (hold : W → Equiv.Perm ι → Prop)
    [∀ w, DecidablePred (hold w)] (c : W → ℝ) (ĉ : W → Equiv.Perm ι → ℝ) {B : ℝ} (hB : 0 ≤ B)
    (hc0 : ∀ w, 0 ≤ c w) (hĉB : ∀ w σ, ĉ w σ ≤ B)
    (hread : ∀ w, ∀ σ ∈ permsOn C, σ t ∉ D w → hold w σ → ĉ w σ ≤ c w) {ν : ℝ}
    (hfail : ∑ w, p w * (((permsOn C).filter fun σ : Equiv.Perm ι =>
        σ t ∉ D w ∧ ¬ hold w σ).card : ℝ) ≤ ν * (permsOn C).card) :
    ∑ w, p w * ((∑ σ ∈ permsOn C, if σ t ∈ D w then 0 else ĉ w σ) / (permsOn C).card) - B * ν ≤
      ∑ w, p w * c w := by
  have hg : (0 : ℝ) < (permsOn C).card := by exact_mod_cast card_permsOn_pos C
  have hw : ∀ w, (∑ σ ∈ permsOn C, if σ t ∈ D w then 0 else ĉ w σ) / (permsOn C).card ≤
      c w + B * ((((permsOn C).filter fun σ : Equiv.Perm ι => σ t ∉ D w ∧ ¬ hold w σ).card : ℝ) /
        (permsOn C).card) := by
    intro w
    have hσ : ∀ σ ∈ permsOn C, (if σ t ∈ D w then 0 else ĉ w σ) ≤
        c w + B * (if σ t ∉ D w ∧ ¬ hold w σ then 1 else 0) := by
      intro σ hσ
      split_ifs with h1 h2 h2
      · linarith [hc0 w]
      · linarith [hc0 w]
      · linarith [hc0 w, hĉB w σ]
      · have hh : hold w σ := by
          by_contra hh
          exact h2 ⟨h1, hh⟩
        linarith [hread w σ hσ h1 hh]
    have hsum := sum_le_sum hσ
    rw [sum_add_distrib, sum_const, nsmul_eq_mul, ← mul_sum, sum_boole] at hsum
    rw [div_le_iff₀ hg, add_mul, mul_assoc, div_mul_cancel₀ _ hg.ne', mul_comm (c w)]
    linarith
  have hsum := sum_le_sum fun w (_ : w ∈ univ) => mul_le_mul_of_nonneg_left (hw w) (hp0 w)
  simp only [mul_add, sum_add_distrib] at hsum
  have hν : ∑ w, p w * ((((permsOn C).filter fun σ : Equiv.Perm ι =>
      σ t ∉ D w ∧ ¬ hold w σ).card : ℝ) / (permsOn C).card) ≤ ν :=
    weighted_share_le C p (fun w σ => σ t ∉ D w ∧ ¬ hold w σ) hfail
  have e2 : ∑ w, p w * (B * ((((permsOn C).filter fun σ : Equiv.Perm ι =>
      σ t ∉ D w ∧ ¬ hold w σ).card : ℝ) / (permsOn C).card)) =
      B * ∑ w, p w * ((((permsOn C).filter fun σ : Equiv.Perm ι =>
        σ t ∉ D w ∧ ¬ hold w σ).card : ℝ) / (permsOn C).card) := by
    rw [mul_sum]
    exact sum_congr rfl fun w _ => by ring
  have e3 := mul_le_mul_of_nonneg_left hν hB
  linarith

/-- **The margin read from the realized ordering.** In ordering `σ` the engine reads a count
`ĉ w σ`, at least `|D ∪ R ∪ M_D|` whenever `σ`'s test unit is kept and `σ`'s bounds hold.
`card_excluded_le_estimated` gives such a reading, `k + |M̂_D| + |D's units the inflated rule
accepts|`, from `σ`'s own estimates (a relabelling of units does not change a count). A left-out
ordering reads nothing. Coverage is at least one, less the orderings' average reading over `|C|`,
less `D`'s share, less `2ν`: `D` is charged inside the count, and again for the orderings that read
nothing. -/
theorem coverage_read_excluding {W Z : Type*} [Fintype W] {C : Finset ι} {t : ι} (ht : t ∈ C)
    {S : (ι → Z) → ι → ℝ} (hS : Equivariant C S) (z : ι → Z) (p : W → ℝ) (hp0 : ∀ w, 0 ≤ p w)
    (hp1 : ∑ w, p w = 1) (ŝ : W → Equiv.Perm ι → ι → ℝ) (ε : W → ι → ℝ) (D : W → Finset ι)
    {ν : ℝ} (hfail : ∑ w, p w * (((permsOn C).filter fun σ : Equiv.Perm ι => σ t ∉ D w ∧
        ¬ ∀ j ∈ C, |ŝ w σ j - S (z ∘ σ) j| ≤ ε w (σ j)).card : ℝ) ≤ ν * (permsOn C).card)
    (k : ℕ) (ĉ : W → Equiv.Perm ι → ℝ) (hĉ0 : ∀ w σ, 0 ≤ ĉ w σ)
    (hread : ∀ w, ∀ σ ∈ permsOn C, σ t ∉ D w →
      (∀ j ∈ C, |ŝ w σ j - S (z ∘ σ) j| ≤ ε w (σ j)) →
      ((C.filter fun u => u ∈ D w ∨ ¬ Accepts k C (S z) u ∨
        Fragile k C (S z) (ε w) u).card : ℝ) ≤ ĉ w σ) :
    1 - ∑ w, p w * ((∑ σ ∈ permsOn C, if σ t ∈ D w then 0 else ĉ w σ / C.card) /
        (permsOn C).card) - ∑ w, p w * (((C.filter fun u => u ∈ D w).card : ℝ) / C.card) -
        2 * ν ≤
      ∑ w, p w * ((((permsOn C).filter fun σ : Equiv.Perm ι =>
        Accepts k C (ŝ w σ) t).card : ℝ) / (permsOn C).card) := by
  have hc : (0 : ℝ) < C.card := by exact_mod_cast card_pos.2 ⟨t, ht⟩
  have h1 := coverage_close_excluding ht hS z p hp0 hp1 ŝ ε D hfail k
  have h2 := le_read_excluding ht p hp0 D
    (fun w σ => ∀ j ∈ C, |ŝ w σ j - S (z ∘ σ) j| ≤ ε w (σ j))
    (fun w => ((C.filter fun u => u ∈ D w ∨ ¬ Accepts k C (S z) u ∨
      Fragile k C (S z) (ε w) u).card : ℝ) / C.card)
    (fun w σ => ĉ w σ / C.card) zero_le_one
    (fun w => div_le_one_of_le₀ (by exact_mod_cast card_filter_le _ _) hc.le)
    (fun w σ => div_nonneg (hĉ0 w σ) hc.le)
    (fun w σ hσ hk hh => div_le_div_of_nonneg_right (hread w σ hσ hk hh) hc.le) hfail
  rw [one_mul] at h2
  linarith

/-- **The equi-tailed set read from the realized ordering**, with a reading at least
`|D ∪ R ∪ M_lo ∪ M_hi|` on kept orderings whose bounds hold
(`card_excluded_le_estimated_twoSided` gives one). -/
theorem coverage_read_twoSided_excluding {W Z : Type*} [Fintype W] {C : Finset ι} {t : ι}
    (ht : t ∈ C) {S : (ι → Z) → ι → ℝ} (hS : Equivariant C S) (z : ι → Z) (p : W → ℝ)
    (hp0 : ∀ w, 0 ≤ p w) (hp1 : ∑ w, p w = 1) (ŝ : W → Equiv.Perm ι → ι → ℝ) (ε : W → ι → ℝ)
    (D : W → Finset ι) {ν : ℝ}
    (hfail : ∑ w, p w * (((permsOn C).filter fun σ : Equiv.Perm ι => σ t ∉ D w ∧
        ¬ ∀ j ∈ C, |ŝ w σ j - S (z ∘ σ) j| ≤ ε w (σ j)).card : ℝ) ≤ ν * (permsOn C).card)
    (k : ℕ) (ĉ : W → Equiv.Perm ι → ℝ) (hĉ0 : ∀ w σ, 0 ≤ ĉ w σ)
    (hread : ∀ w, ∀ σ ∈ permsOn C, σ t ∉ D w →
      (∀ j ∈ C, |ŝ w σ j - S (z ∘ σ) j| ≤ ε w (σ j)) →
      ((C.filter fun u => u ∈ D w ∨ ¬ AcceptsTwoSided k C (S z) u ∨
        FragileTwoSided k C (S z) (ε w) u).card : ℝ) ≤ ĉ w σ) :
    1 - ∑ w, p w * ((∑ σ ∈ permsOn C, if σ t ∈ D w then 0 else ĉ w σ / C.card) /
        (permsOn C).card) - ∑ w, p w * (((C.filter fun u => u ∈ D w).card : ℝ) / C.card) -
        2 * ν ≤
      ∑ w, p w * ((((permsOn C).filter fun σ : Equiv.Perm ι =>
        AcceptsTwoSided k C (ŝ w σ) t).card : ℝ) / (permsOn C).card) := by
  have hc : (0 : ℝ) < C.card := by exact_mod_cast card_pos.2 ⟨t, ht⟩
  have h1 := coverage_close_twoSided_excluding ht hS z p hp0 hp1 ŝ ε D hfail k
  have h2 := le_read_excluding ht p hp0 D
    (fun w σ => ∀ j ∈ C, |ŝ w σ j - S (z ∘ σ) j| ≤ ε w (σ j))
    (fun w => ((C.filter fun u => u ∈ D w ∨ ¬ AcceptsTwoSided k C (S z) u ∨
      FragileTwoSided k C (S z) (ε w) u).card : ℝ) / C.card)
    (fun w σ => ĉ w σ / C.card) zero_le_one
    (fun w => div_le_one_of_le₀ (by exact_mod_cast card_filter_le _ _) hc.le)
    (fun w σ => div_nonneg (hĉ0 w σ) hc.le)
    (fun w σ hσ hk hh => div_le_div_of_nonneg_right (hread w σ hσ hk hh) hc.le) hfail
  rw [one_mul] at h2
  linarith

/-- **The PIT read from the realized ordering.** Readings `up` (from above) and `lo` (from below) of
`pit_excluding`'s counts on kept orderings whose bounds hold (`card_pit_excluding_estimated` gives
both), nothing on a left-out ordering. `P(PIT ≤ k/|C|)` is at least the average lower reading
less `2ν`, and at most the average upper reading plus `D`'s share plus `2ν`. -/
theorem pit_read_excluding {W Z : Type*} [Fintype W] {C : Finset ι} {t : ι} (ht : t ∈ C)
    {S : (ι → Z) → ι → ℝ} (hS : Equivariant C S) (z : ι → Z) (p : W → ℝ) (hp0 : ∀ w, 0 ≤ p w)
    (ŝ : W → Equiv.Perm ι → ι → ℝ) (ε : W → ι → ℝ) (D : W → Finset ι) {ν : ℝ}
    (hfail : ∑ w, p w * (((permsOn C).filter fun σ : Equiv.Perm ι => σ t ∉ D w ∧
        ¬ ∀ j ∈ C, |ŝ w σ j - S (z ∘ σ) j| ≤ ε w (σ j)).card : ℝ) ≤ ν * (permsOn C).card)
    (k : ℕ) (up lo : W → Equiv.Perm ι → ℝ) (hup0 : ∀ w σ, 0 ≤ up w σ)
    (hlo1 : ∀ w σ, lo w σ ≤ C.card)
    (hreadA : ∀ w, ∀ σ ∈ permsOn C, σ t ∉ D w →
      (∀ j ∈ C, |ŝ w σ j - S (z ∘ σ) j| ≤ ε w (σ j)) →
      ((C.filter fun u => u ∈ D w ∨ shiftedBelow C (S z) (ε w) u < k).card : ℝ) ≤ up w σ)
    (hreadB : ∀ w, ∀ σ ∈ permsOn C, σ t ∉ D w →
      (∀ j ∈ C, |ŝ w σ j - S (z ∘ σ) j| ≤ ε w (σ j)) →
      lo w σ ≤ ((C.filter fun u => u ∉ D w ∧ shiftedAtMost C (S z) (ε w) u ≤ k).card : ℝ)) :
    ∑ w, p w * ((∑ σ ∈ permsOn C, if σ t ∈ D w then 0 else lo w σ / C.card) /
        (permsOn C).card) - 2 * ν ≤
      ∑ w, p w * ((∑ σ ∈ permsOn C, pitRankLe C (ŝ w σ) t k) / (permsOn C).card) ∧
    ∑ w, p w * ((∑ σ ∈ permsOn C, pitRankLe C (ŝ w σ) t k) / (permsOn C).card) ≤
      ∑ w, p w * ((∑ σ ∈ permsOn C, if σ t ∈ D w then 0 else up w σ / C.card) /
        (permsOn C).card) + ∑ w, p w * (((C.filter fun u => u ∈ D w).card : ℝ) / C.card) +
        2 * ν := by
  have hc : (0 : ℝ) < C.card := by exact_mod_cast card_pos.2 ⟨t, ht⟩
  obtain ⟨hlo, hup⟩ := pit_excluding ht hS z p hp0 ŝ ε D hfail k
  constructor
  · have h2 := read_le_excluding p hp0 D
      (fun w σ => ∀ j ∈ C, |ŝ w σ j - S (z ∘ σ) j| ≤ ε w (σ j))
      (fun w => ((C.filter fun u => u ∉ D w ∧ shiftedAtMost C (S z) (ε w) u ≤ k).card : ℝ) /
        C.card)
      (fun w σ => lo w σ / C.card) zero_le_one (fun w => div_nonneg (Nat.cast_nonneg _) hc.le)
      (fun w σ => div_le_one_of_le₀ (hlo1 w σ) hc.le)
      (fun w σ hσ hk hh => div_le_div_of_nonneg_right (hreadB w σ hσ hk hh) hc.le) hfail
    rw [one_mul] at h2
    linarith
  · have h2 := le_read_excluding ht p hp0 D
      (fun w σ => ∀ j ∈ C, |ŝ w σ j - S (z ∘ σ) j| ≤ ε w (σ j))
      (fun w => ((C.filter fun u => u ∈ D w ∨ shiftedBelow C (S z) (ε w) u < k).card : ℝ) /
        C.card)
      (fun w σ => up w σ / C.card) zero_le_one
      (fun w => div_le_one_of_le₀ (by exact_mod_cast card_filter_le _ _) hc.le)
      (fun w σ => div_nonneg (hup0 w σ) hc.le)
      (fun w σ hσ hk hh => div_le_div_of_nonneg_right (hreadA w σ hσ hk hh) hc.le) hfail
    rw [one_mul] at h2
    linarith

end ReadingExcluded

/-! ### Slice 10's grid reading: the PIT between grid values, and the margin over the grid -/

section GridReading

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- The chance that `V`, uniform on `(0, 1]`, is at most `a`. -/
noncomputable def unitClamp (a : ℝ) : ℝ := max 0 (min 1 a)

omit [Fintype ι] [DecidableEq ι] in
/-- **The full step's PIT CDF**: `P(PIT ≤ x)` given the scores, the PIT being `(R − 1 + V)/|C|`,
with `R` the randomized rank (`tieRankProb`) and `V` uniform on `(0, 1]`. -/
noncomputable def pitCdf (C : Finset ι) (s : ι → ℝ) (u : ι) (x : ℝ) : ℝ :=
  ∑ r ∈ Finset.Icc 1 C.card, tieRankProb C s u r * unitClamp (C.card * x - r + 1)

omit [Fintype ι] [DecidableEq ι] in
/-- `P(R ≤ k)` as a sum over the ranks `1..|C|`. -/
theorem pitRankLe_eq_sum_ite (C : Finset ι) (s : ι → ℝ) (u : ι) (k : ℕ) :
    pitRankLe C s u k =
      ∑ r ∈ Finset.Icc 1 C.card, tieRankProb C s u r * (if r ≤ k then 1 else 0) := by
  classical
  have hz : ∀ r, C.card < r → tieRankProb C s u r = 0 := fun r hr => by
    unfold tieRankProb
    rw [if_neg]
    rintro ⟨-, h⟩
    have := card_filter_le C (fun j => s j ≤ s u)
    unfold rankLe at h
    omega
  unfold pitRankLe
  simp only [mul_ite, mul_one, mul_zero]
  rw [← sum_filter]
  refine (sum_subset (fun r hr => ?_) fun r hr hn => hz r ?_).symm
  · rw [mem_filter, Finset.mem_Icc] at hr
    exact Finset.mem_Icc.2 ⟨hr.1.1, hr.2⟩
  · rw [Finset.mem_Icc] at hr
    by_contra hle
    exact hn (mem_filter.2 ⟨Finset.mem_Icc.2 ⟨hr.1, not_lt.1 hle⟩, hr.2⟩)

omit [Fintype ι] [DecidableEq ι] in
/-- **The PIT is linear between grid values**: at `x = (k + θ)/|C|`, `θ ∈ [0, 1]`,
`P(PIT ≤ x) = (1 − θ) P(R ≤ k) + θ P(R ≤ k + 1)`. -/
theorem pitCdf_interp (C : Finset ι) (s : ι → ℝ) {u : ι} (hu : u ∈ C) (k : ℕ) {θ : ℝ}
    (h0 : 0 ≤ θ) (h1 : θ ≤ 1) :
    pitCdf C s u ((k + θ) / C.card) =
      (1 - θ) * pitRankLe C s u k + θ * pitRankLe C s u (k + 1) := by
  have hc : (0 : ℝ) < C.card := by exact_mod_cast card_pos.2 ⟨u, hu⟩
  have e : (C.card : ℝ) * ((k + θ) / C.card) = k + θ := by field_simp
  rw [pitRankLe_eq_sum_ite, pitRankLe_eq_sum_ite, mul_sum, mul_sum, ← sum_add_distrib]
  unfold pitCdf
  rw [e]
  refine sum_congr rfl fun r _ => ?_
  rcases lt_trichotomy r (k + 1) with hlt | heq | hgt
  · have hr : r ≤ k := by omega
    have hrr : (r : ℝ) ≤ k := by exact_mod_cast hr
    have hcl : unitClamp ((k : ℝ) + θ - r + 1) = 1 := by
      unfold unitClamp
      rw [min_eq_left (by linarith), max_eq_right zero_le_one]
    rw [hcl, if_pos hr, if_pos (by omega)]
    ring
  · subst heq
    have hcl : unitClamp ((k : ℝ) + θ - ((k + 1 : ℕ) : ℝ) + 1) = θ := by
      unfold unitClamp
      push_cast
      rw [show (k : ℝ) + θ - (k + 1) + 1 = θ by ring, min_eq_right h1, max_eq_right h0]
    rw [hcl, if_neg (by omega), if_pos le_rfl]
    ring
  · have hrr : (k : ℝ) + 2 ≤ r := by exact_mod_cast (show k + 2 ≤ r by omega)
    have hcl : unitClamp ((k : ℝ) + θ - r + 1) = 0 := by
      unfold unitClamp
      rw [min_eq_right (by linarith), max_eq_left (by linarith)]
    rw [hcl, if_neg (by omega), if_neg (by omega)]
    ring

/-- The interpolation's departure from the diagonal is at most the larger at its ends. -/
theorem interp_dev_le (k : ℕ) {m θ A B : ℝ} (hm : 0 < m) (h0 : 0 ≤ θ) (h1 : θ ≤ 1) :
    |(1 - θ) * A + θ * B - (k + θ) / m| ≤ max |A - k / m| |B - (k + 1) / m| := by
  have e : (1 - θ) * A + θ * B - (k + θ) / m =
      (1 - θ) * (A - k / m) + θ * (B - (k + 1) / m) := by
    field_simp
    ring
  rw [e]
  calc |(1 - θ) * (A - k / m) + θ * (B - (k + 1) / m)|
      ≤ |(1 - θ) * (A - k / m)| + |θ * (B - (k + 1) / m)| := abs_add_le _ _
    _ = (1 - θ) * |A - k / m| + θ * |B - (k + 1) / m| := by
        rw [abs_mul, abs_mul, abs_of_nonneg (by linarith : (0 : ℝ) ≤ 1 - θ), abs_of_nonneg h0]
    _ ≤ (1 - θ) * max |A - k / m| |B - (k + 1) / m| +
          θ * max |A - k / m| |B - (k + 1) / m| :=
        add_le_add (mul_le_mul_of_nonneg_left (le_max_left _ _) (by linarith))
          (mul_le_mul_of_nonneg_left (le_max_right _ _) h0)
    _ = max |A - k / m| |B - (k + 1) / m| := by ring

/-- **The PIT averaged over the draws and the orderings is linear between grid values too**, so
`sup_x |P(PIT ≤ x) − x|` is the largest departure at the grid values `k/|C|`
(`interp_dev_le`), and bounds on the grid values (`pit_excluding`, `pit_read_excluding`)
interpolate. -/
theorem avg_pitCdf_interp {W : Type*} [Fintype W] {C : Finset ι} {t : ι} (ht : t ∈ C)
    (p : W → ℝ) (ŝ : W → Equiv.Perm ι → ι → ℝ) (k : ℕ) {θ : ℝ} (h0 : 0 ≤ θ) (h1 : θ ≤ 1) :
    ∑ w, p w * ((∑ σ ∈ permsOn C, pitCdf C (ŝ w σ) t ((k + θ) / C.card)) / (permsOn C).card) =
      (1 - θ) * ∑ w, p w * ((∑ σ ∈ permsOn C, pitRankLe C (ŝ w σ) t k) / (permsOn C).card) +
        θ * ∑ w, p w * ((∑ σ ∈ permsOn C, pitRankLe C (ŝ w σ) t (k + 1)) /
          (permsOn C).card) := by
  simp only [pitCdf_interp _ _ ht k h0 h1]
  rw [mul_sum, mul_sum, ← sum_add_distrib]
  refine sum_congr rfl fun w _ => ?_
  rw [sum_add_distrib, ← mul_sum, ← mul_sum]
  ring

/-- The largest of `f` over the arrangements that put each grid value at the test position. -/
noncomputable def gridMax {Z : Type*} (G : Finset Z) (hG : G.Nonempty) (f : (ι → Z) → ℝ)
    (t : ι) (x : ι → Z) : ℝ :=
  G.sup' hG fun y => f (Function.update x t y)

omit [Fintype ι] in
/-- **The reading at the true value is one of the grid's**, when the grid holds that value. -/
theorem le_gridMax {Z : Type*} {G : Finset Z} (hG : G.Nonempty) (f : (ι → Z) → ℝ) {t : ι}
    {x : ι → Z} (hx : x t ∈ G) : f x ≤ gridMax G hG f t x := by
  unfold gridMax
  have h := le_sup' (fun y => f (Function.update x t y)) hx
  simp only [Function.update_eq_self] at h
  exact h

omit [Fintype ι] in
/-- **The grid's report is free of the test unit's own value.** -/
theorem gridMax_update {Z : Type*} {G : Finset Z} (hG : G.Nonempty) (f : (ι → Z) → ℝ) (t : ι)
    (x : ι → Z) (y : Z) : gridMax G hG f t (Function.update x t y) = gridMax G hG f t x := by
  unfold gridMax
  simp only [Function.update_idem]

omit [Fintype ι] in
/-- The engine's report over a grid: the largest reading at a value whose turn is kept, and the
largest flagged count, both over the grid's values. -/
noncomputable def gridReport {Z : Type*} (C : Finset ι) (t : ι) (grid : (ι → Z) → Finset Z)
    (hne : ∀ x, (grid x).Nonempty) (eng : (ι → Z) → ℝ) (flagged : (ι → Z) → Finset ι)
    (x : ι → Z) : ℝ :=
  gridMax (grid x) (hne x) (fun x' => if t ∈ flagged x' then 0 else eng x') t x +
    gridMax (grid x) (hne x) (fun x' => ((C.filter fun u => u ∈ flagged x').card : ℝ)) t x

/-- **The margin over the grid.** The engine reads `eng w x` from an arrangement `x` (the data with
the test unit's value at `t`), and flags the positions `flagged w x`. The realized ordering's flag
on its test unit is the data set's `D`, and the flagged count covers `|D|`. The grid holds the
value the test unit takes, possibly built from the rest of the data (the test unit's covariate).
Then coverage is at least one, less the average over the orderings of the grid's report
(`gridReport`: the largest kept reading plus the largest flagged count, over the grid's values) over
`|C|`, less `2ν`. With the grid fixed, the report is free of the test unit's own value
(`gridMax_update`), so it is read from the data. A grid that misses the value the test unit takes
can miss its reading. -/
theorem coverage_grid_excluding {W Z : Type*} [Fintype W] {C : Finset ι} {t : ι} (ht : t ∈ C)
    {S : (ι → Z) → ι → ℝ} (hS : Equivariant C S) (z : ι → Z) (p : W → ℝ) (hp0 : ∀ w, 0 ≤ p w)
    (hp1 : ∑ w, p w = 1) (ŝ : W → Equiv.Perm ι → ι → ℝ) (ε : W → ι → ℝ) (D : W → Finset ι)
    {ν : ℝ} (hfail : ∑ w, p w * (((permsOn C).filter fun σ : Equiv.Perm ι => σ t ∉ D w ∧
        ¬ ∀ j ∈ C, |ŝ w σ j - S (z ∘ σ) j| ≤ ε w (σ j)).card : ℝ) ≤ ν * (permsOn C).card)
    (k : ℕ) (grid : (ι → Z) → Finset Z) (hne : ∀ x, (grid x).Nonempty)
    (hgrid : ∀ σ ∈ permsOn C, (z ∘ σ) t ∈ grid (z ∘ σ)) (eng : W → (ι → Z) → ℝ)
    (heng0 : ∀ w x, 0 ≤ eng w x) (flagged : W → (ι → Z) → Finset ι)
    (hflag : ∀ w, ∀ σ ∈ permsOn C, (σ t ∈ D w ↔ t ∈ flagged w (z ∘ σ)))
    (hcount : ∀ w, ∀ σ ∈ permsOn C, (C.filter fun u => u ∈ D w).card ≤
      (C.filter fun u => u ∈ flagged w (z ∘ σ)).card)
    (hread : ∀ w, ∀ σ ∈ permsOn C, t ∉ flagged w (z ∘ σ) →
      (∀ j ∈ C, |ŝ w σ j - S (z ∘ σ) j| ≤ ε w (σ j)) →
      ((C.filter fun u => u ∈ D w ∨ ¬ Accepts k C (S z) u ∨
        Fragile k C (S z) (ε w) u).card : ℝ) ≤ eng w (z ∘ σ)) :
    1 - ∑ w, p w * ((∑ σ ∈ permsOn C,
        gridReport C t grid hne (eng w) (flagged w) (z ∘ σ) / C.card) / (permsOn C).card) -
        2 * ν ≤
      ∑ w, p w * ((((permsOn C).filter fun σ : Equiv.Perm ι =>
        Accepts k C (ŝ w σ) t).card : ℝ) / (permsOn C).card) := by
  have hc : (0 : ℝ) < C.card := by exact_mod_cast card_pos.2 ⟨t, ht⟩
  have hg : (0 : ℝ) < (permsOn C).card := by exact_mod_cast card_permsOn_pos C
  have key := coverage_read_excluding ht hS z p hp0 hp1 ŝ ε D hfail k
    (fun w σ => eng w (z ∘ σ)) (fun w σ => heng0 w _)
    (fun w σ hσ hk hh => hread w σ hσ (fun hf => hk ((hflag w σ hσ).2 hf)) hh)
  have hb : ∀ w, (∑ σ ∈ permsOn C, if σ t ∈ D w then 0 else eng w (z ∘ σ) / C.card) /
      (permsOn C).card + ((C.filter fun u => u ∈ D w).card : ℝ) / C.card ≤
      (∑ σ ∈ permsOn C, gridReport C t grid hne (eng w) (flagged w) (z ∘ σ) / C.card) /
        (permsOn C).card := by
    intro w
    have e : ((C.filter fun u => u ∈ D w).card : ℝ) / C.card =
        (∑ _σ ∈ permsOn C, ((C.filter fun u => u ∈ D w).card : ℝ) / C.card) /
          (permsOn C).card := by
      rw [sum_const, nsmul_eq_mul, mul_div_cancel_left₀ _ hg.ne']
    rw [e, ← add_div, ← sum_add_distrib]
    refine div_le_div_of_nonneg_right (sum_le_sum fun σ hσ => ?_) hg.le
    have h1 := le_gridMax (hne (z ∘ σ)) (fun x' => if t ∈ flagged w x' then 0 else eng w x')
      (hgrid σ hσ)
    have h2 := le_gridMax (hne (z ∘ σ))
      (fun x' => ((C.filter fun u => u ∈ flagged w x').card : ℝ)) (hgrid σ hσ)
    have h3 : ((C.filter fun u => u ∈ D w).card : ℝ) ≤
        ((C.filter fun u => u ∈ flagged w (z ∘ σ)).card : ℝ) := by
      exact_mod_cast hcount w σ hσ
    unfold gridReport
    rw [add_div]
    refine add_le_add ?_ (div_le_div_of_nonneg_right (h3.trans h2) hc.le)
    by_cases hD : σ t ∈ D w
    · have hf := (hflag w σ hσ).1 hD
      simp only [hf, ↓reduceIte] at h1
      rw [if_pos hD]
      exact div_nonneg h1 hc.le
    · have hf : t ∉ flagged w (z ∘ σ) := fun h => hD ((hflag w σ hσ).2 h)
      simp only [hf, ↓reduceIte] at h1
      rw [if_neg hD]
      exact div_le_div_of_nonneg_right h1 hc.le
  have hsum := sum_le_sum fun w (_ : w ∈ univ) => mul_le_mul_of_nonneg_left (hb w) (hp0 w)
  simp only [mul_add, sum_add_distrib] at hsum
  linarith

end GridReading

/-! ### Slice 13: the tournament correction (Bhattacharyya, Zhang and Barber, arXiv 2605.29200) -/

section Tournament

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- **The tournament correction's pair fit.** The fit to the units `T`, then the scored value with
the units `T` leaves out handed over as an unordered multiset of values. On the leave-two-out set
`univ \ {i, j}`, `pairScore` of it is their `s_approx(z_i; D_{−ij} + {z_i, z_j})`, a conformity in
place of their nonconformity. -/
def tournFit {Z β α : Type*} (fit : Finset ι → (ι → Z) → β) (g : β → Z → Multiset Z → α)
    (T : Finset ι) (z : ι → Z) (p : Z) : α :=
  g (fit T z) p ((univ \ T).val.map z)

theorem compl_image_val_map {Z : Type*} (σ : Equiv.Perm ι) (T : Finset ι) (z : ι → Z) :
    (univ \ T.image σ).val.map z = (univ \ T).val.map (z ∘ σ) := by
  have h : univ \ T.image σ = (univ \ T).map σ.toEmbedding := by
    ext x
    simp only [mem_sdiff, mem_univ, true_and, mem_image, mem_map, Equiv.coe_toEmbedding]
    constructor
    · intro hx
      exact ⟨σ.symm x, fun h => hx ⟨σ.symm x, h, σ.apply_symm_apply x⟩, σ.apply_symm_apply x⟩
    · rintro ⟨a, ha, rfl⟩ ⟨b, hb, hba⟩
      exact ha (σ.injective hba ▸ hb)
  rw [h, Finset.map_val, Multiset.map_map]
  rfl

/-- The pair fit inherits the fit's relabeling: the pair is passed as values, not positions. -/
theorem tournFit_comp {Z β α : Type*} {fit fit' : Finset ι → (ι → Z) → β}
    (g : β → Z → Multiset Z → α) {σ : Equiv.Perm ι} {z : ι → Z} {T : Finset ι}
    (h : fit' T (z ∘ σ) = fit (T.image σ) z) :
    tournFit fit' g T (z ∘ σ) = tournFit fit g (T.image σ) z := by
  funext p
  simp only [tournFit, h, compl_image_val_map]

/-- **The tournament correction covers at least `1 − (2k − 1)/|C|`**, so `1 − 2α` at
`k = ⌊α|C|⌋` (`one_sub_two_mul_le`): their Theorem 1, randomized scores included (their Theorem
3). Each pair's fit to the rest is symmetric in its training units given its randomness, up to a
relabeling `act` (seeds that travel with the units qualify); the pair enters only through `g`, as
an unordered multiset of values. -/
theorem tournament_coverage {W Z β α : Type*} [Fintype W] [Fintype Z] [LinearOrder α]
    {C : Finset ι} {t : ι} (ht : t ∈ C) (P : W → (ι → Z) → ℝ) (hP0 : ∀ w z, 0 ≤ P w z)
    (hP1 : ∑ w, ∑ z, P w z = 1) (act : Equiv.Perm ι → W ≃ W)
    (hP : ∀ σ ∈ permsOn C, ∀ w z, P (act σ w) (z ∘ σ) = P w z)
    (fit : W → Finset ι → (ι → Z) → β) (g : β → Z → Multiset Z → α)
    (hfit : ∀ σ ∈ permsOn C, ∀ w T z, fit (act σ w) T (z ∘ σ) = fit w (T.image σ) z) (k : ℕ) :
    1 - ((2 * k - 1 : ℕ) : ℝ) / C.card ≤
      ∑ w, ∑ z, P w z *
        if AcceptsPairwise k C (pairScore (tournFit (fit w) g) z) t then 1 else 0 :=
  crossFit_coverage ht P hP0 hP1 act hP (fun w => tournFit (fit w) g)
    (fun σ hσ w T z => tournFit_comp g (hfit σ hσ w T z)) k

/-- **Each pair's chain need not be symmetric, only read its units' order.** -/
theorem tournFit_orderReading {n : ℕ} {Z β α : Type*}
    {fit : Finset (Fin (n + 1)) → (Fin (n + 1) → Z) → β} (g : β → Z → Multiset Z → α)
    (h : OrderReading fit) : OrderReading (tournFit fit g) :=
  fun σ T z hmono => tournFit_comp g (h σ T z hmono)

/-- **The tournament correction with fits that read their units' order**: a chain sweeping its
training units in index order with a fixed seed, the left-out places closed up, need be neither
symmetric nor mixed. The test unit is unit `0`. -/
theorem tournament_orderReading {n : ℕ} {W Z β α : Type*} [Fintype W] [Fintype Z] [LinearOrder α]
    (P : W → (Fin (n + 1) → Z) → ℝ) (hP0 : ∀ w z, 0 ≤ P w z) (hP1 : ∑ w, ∑ z, P w z = 1)
    (hP : ∀ σ : Equiv.Perm (Fin (n + 1)), ∀ w z, P w (z ∘ σ) = P w z)
    (fit : W → Finset (Fin (n + 1)) → (Fin (n + 1) → Z) → β) (g : β → Z → Multiset Z → α)
    (hfit : ∀ w, OrderReading (fit w)) (k : ℕ) :
    1 - ((2 * k - 1 : ℕ) : ℝ) / (n + 1) ≤
      ∑ w, ∑ z, P w z *
        if AcceptsPairwise k univ (pairScore (tournFit (fit w) g) z) 0 then 1 else 0 :=
  coverage_orderReading P hP0 hP1 hP (fun w => tournFit (fit w) g)
    (fun w => tournFit_orderReading g (hfit w)) k

omit [Fintype ι] [DecidableEq ι] in
/-- **Their threshold is the cross-fitted rule's level.** They refuse when the test unit wins at
least `(1 − α)m` games, a win being a strictly higher nonconformity: a strictly lower conformity, so
ties win nothing, as in `strangerThan`. With `c` wins that is `c + ⌊αm⌋ ≥ m`, the cross-fitted rule
at `k = ⌊αm⌋`. -/
theorem tournament_count_iff {c m : ℕ} {α : ℝ} (hα : 0 ≤ α) :
    (c : ℝ) < (1 - α) * m ↔ c + ⌊α * m⌋₊ < m := by
  have h0 : (0 : ℝ) ≤ α * m := mul_nonneg hα (Nat.cast_nonneg _)
  constructor
  · intro h
    have hcm : c < m := by
      by_contra hc
      have : (m : ℝ) ≤ c := by exact_mod_cast not_lt.1 hc
      nlinarith
    have h2 : α * m < ((m - c : ℕ) : ℝ) := by
      rw [Nat.cast_sub hcm.le]
      linarith
    have := (Nat.floor_lt h0).2 h2
    omega
  · intro h
    have h1 : ⌊α * m⌋₊ < m - c := by omega
    have h2 := (Nat.floor_lt h0).1 h1
    rw [Nat.cast_sub (by omega)] at h2
    linarith

/-! #### The `2ε`-inflated rule under stability -/

omit [Fintype ι] in
/-- The units of `C` that `u` is stranger than by more than `δ`: `s u j + δ < s j u`. -/
noncomputable def strangerThanBy (C : Finset ι) (s : ι → ι → ℝ) (δ : ℝ) (u : ι) : Finset ι :=
  C.filter fun j => s u j + δ < s j u

omit [Fintype ι] in
/-- **The inflated tournament rule**: their `2ε`-inflated set at `δ = 2ε`. -/
def AcceptsPairwiseBy (k : ℕ) (C : Finset ι) (s : ι → ι → ℝ) (δ : ℝ) (u : ι) : Prop :=
  (strangerThanBy C s δ u).card + k < C.card

omit [Fintype ι] in
noncomputable instance (k : ℕ) (C : Finset ι) (s : ι → ι → ℝ) (δ : ℝ) (u : ι) :
    Decidable (AcceptsPairwiseBy k C s δ u) :=
  inferInstanceAs (Decidable (_ < _))

omit [Fintype ι] in
/-- The ordered pairs whose score strays from a reference `R` by more than `ε`. -/
noncomputable def strayCount (C : Finset ι) (s : ι → ι → ℝ) (R : ι → ℝ) (ε : ℝ) : ℕ :=
  ∑ u ∈ C, (C.filter fun j => j ≠ u ∧ ε < |s u j - R u|).card

omit [Fintype ι] in
/-- The games of `u` in which either side strays. -/
noncomputable def pairFails (C : Finset ι) (s : ι → ι → ℝ) (R : ι → ℝ) (ε : ℝ) (u : ι) : ℕ :=
  (C.filter fun j => j ≠ u ∧ (ε < |s u j - R u| ∨ ε < |s j u - R j|)).card

omit [Fintype ι] in
/-- **The inflated rule's refusals, counted.** For any reference scores `R` and every `a`, at most
`k + a + 2F/(a + 1)` units are refused, `F` the stray pairs. Where neither side of a game strays,
an inflated win needs `R u < R j`; so a refused unit either sits among the `k + a` lowest by `R`, or
has more than `a` stray games, and Markov counts those. -/
theorem card_refused_by_le (C : Finset ι) (s : ι → ι → ℝ) (R : ι → ℝ) {ε : ℝ} (hε : 0 ≤ ε)
    (k a : ℕ) :
    ((C.filter fun u => ¬ AcceptsPairwiseBy k C s (2 * ε) u).card : ℝ) ≤
      (k + a : ℝ) + 2 * (strayCount C s R ε : ℝ) / (a + 1) := by
  classical
  have hunit : ∀ u ∈ C, ¬ AcceptsPairwiseBy k C s (2 * ε) u →
      rankLe C R u ≤ k + pairFails C s R ε u := by
    intro u _ hna
    rw [AcceptsPairwiseBy, not_lt] at hna
    have h1 : (strangerThanBy C s (2 * ε) u).card ≤
        (C.filter fun j => R u < R j).card + pairFails C s R ε u := by
      refine (card_le_card fun j hj => ?_).trans (card_union_le _ _)
      rw [strangerThanBy, mem_filter] at hj
      rw [mem_union, mem_filter, mem_filter]
      by_cases hf : j ≠ u ∧ (ε < |s u j - R u| ∨ ε < |s j u - R j|)
      · exact Or.inr ⟨hj.1, hf⟩
      · left
        refine ⟨hj.1, ?_⟩
        have hju : j ≠ u := by
          rintro rfl
          linarith [hj.2]
        have hnot : ¬ (ε < |s u j - R u| ∨ ε < |s j u - R j|) := fun h => hf ⟨hju, h⟩
        rw [not_or, not_lt, not_lt] at hnot
        have h1 := abs_le.1 hnot.1
        have h2 := abs_le.1 hnot.2
        linarith [hj.2]
    have h2 : (C.filter fun j => R u < R j).card + rankLe C R u = C.card := by
      unfold rankLe
      have := C.card_filter_add_card_filter_not (fun j => R j ≤ R u)
      rw [filter_congr fun j _ => show ¬ R j ≤ R u ↔ R u < R j from not_le] at this
      omega
    omega
  have hcover : (C.filter fun u => ¬ AcceptsPairwiseBy k C s (2 * ε) u) ⊆
      (C.filter fun u => rankLe C R u ≤ k + a) ∪ C.filter fun u => a + 1 ≤ pairFails C s R ε u := by
    intro u hu
    rw [mem_filter] at hu
    rw [mem_union, mem_filter, mem_filter]
    have := hunit u hu.1 hu.2
    by_cases hfa : pairFails C s R ε u ≤ a
    · exact Or.inl ⟨hu.1, by omega⟩
    · exact Or.inr ⟨hu.1, by omega⟩
  have hmarkov : (a + 1) * (C.filter fun u => a + 1 ≤ pairFails C s R ε u).card ≤
      ∑ u ∈ C, pairFails C s R ε u := by
    calc (a + 1) * (C.filter fun u => a + 1 ≤ pairFails C s R ε u).card
        = ∑ _u ∈ C.filter (fun u => a + 1 ≤ pairFails C s R ε u), (a + 1) := by
          rw [sum_const, smul_eq_mul, mul_comm]
      _ ≤ ∑ u ∈ C.filter (fun u => a + 1 ≤ pairFails C s R ε u), pairFails C s R ε u :=
          sum_le_sum fun u hu => (mem_filter.1 hu).2
      _ ≤ ∑ u ∈ C, pairFails C s R ε u :=
          sum_le_sum_of_subset_of_nonneg (filter_subset _ _) fun _ _ _ => Nat.zero_le _
  have hfails : ∀ u ∈ C, pairFails C s R ε u ≤
      (C.filter fun j => j ≠ u ∧ ε < |s u j - R u|).card +
        (C.filter fun j => j ≠ u ∧ ε < |s j u - R j|).card := by
    intro u _
    refine (card_le_card fun j hj => ?_).trans (card_union_le _ _)
    rw [mem_filter] at hj
    rw [mem_union, mem_filter, mem_filter]
    rcases hj.2.2 with h | h
    · exact Or.inl ⟨hj.1, hj.2.1, h⟩
    · exact Or.inr ⟨hj.1, hj.2.1, h⟩
  have hB : ∑ u ∈ C, (C.filter fun j => j ≠ u ∧ ε < |s j u - R j|).card = strayCount C s R ε := by
    unfold strayCount
    simp only [card_filter]
    rw [sum_comm]
    refine sum_congr rfl fun u _ => sum_congr rfl fun j _ => ?_
    by_cases h : u = j
    · simp [h]
    · simp [h, Ne.symm h]
  have hsumF : ∑ u ∈ C, pairFails C s R ε u ≤ 2 * strayCount C s R ε := by
    calc ∑ u ∈ C, pairFails C s R ε u
        ≤ ∑ u ∈ C, ((C.filter fun j => j ≠ u ∧ ε < |s u j - R u|).card +
            (C.filter fun j => j ≠ u ∧ ε < |s j u - R j|).card) := sum_le_sum hfails
      _ = strayCount C s R ε + strayCount C s R ε := by
          rw [sum_add_distrib, hB]
          rfl
      _ = 2 * strayCount C s R ε := by ring
  have htot := (card_le_card hcover).trans (card_union_le _ _)
  have hc1 := card_rankLe_le C R (k + a)
  have hN : ((C.filter fun u => a + 1 ≤ pairFails C s R ε u).card : ℝ) ≤
      2 * (strayCount C s R ε : ℝ) / (a + 1) := by
    rw [le_div_iff₀ (by positivity : (0 : ℝ) < a + 1)]
    have : (C.filter fun u => a + 1 ≤ pairFails C s R ε u).card * (a + 1) ≤
        2 * strayCount C s R ε := by
      rw [mul_comm]
      exact hmarkov.trans hsumF
    exact_mod_cast this
  have h1 : ((C.filter fun u => ¬ AcceptsPairwiseBy k C s (2 * ε) u).card : ℝ) ≤
      ((C.filter fun u => rankLe C R u ≤ k + a).card : ℝ) +
        ((C.filter fun u => a + 1 ≤ pairFails C s R ε u).card : ℝ) := by
    exact_mod_cast htot
  have h2 : ((C.filter fun u => rankLe C R u ≤ k + a).card : ℝ) ≤ k + a := by
    exact_mod_cast hc1
  linarith

theorem acceptsPairwiseBy_comp {k : ℕ} {C : Finset ι} {σ : Equiv.Perm ι} (hσ : σ ∈ permsOn C)
    (s : ι → ι → ℝ) (δ : ℝ) (u : ι) :
    AcceptsPairwiseBy k C (fun i j => s (σ i) (σ j)) δ u ↔ AcceptsPairwiseBy k C s δ (σ u) := by
  have h : (strangerThanBy C (fun i j => s (σ i) (σ j)) δ u).card =
      (strangerThanBy C s δ (σ u)).card := by
    apply card_equiv σ
    intro j
    simp only [strangerThanBy, mem_filter, mem_iff_of_mem_permsOn hσ]
  simp only [AcceptsPairwiseBy, h]

/-- **The `2ε`-inflated tournament rule under stability**, per data set (their Theorem 2). Over the
orderings of a data set, with a pairwise-equivariant score, the test unit is accepted at least
`1 − (k + a + 2F/(a + 1))/|C|` of the time, for every `a` and **any** reference scores `R` of the
data set, `F` the ordered pairs whose score strays from `R` by more than `ε`. Under an exchangeable
law in which each ordered pair strays with chance at most `ν` (their assumption, with their
one-out scores for `R`), `E F ≤ |C|(|C| − 1)ν`, so the miscoverage is at most
`k/|C| + a/|C| + 2(|C| − 1)ν/(a + 1)` for every `a`; optimized over `a` this is about
`α + 2√(2ν)`, below their `α + 3√ν`. -/
theorem coverage_tournament_stable {Z : Type*} {C : Finset ι} {t : ι} (ht : t ∈ C)
    {S : (ι → Z) → ι → ι → ℝ} (hS : PairEquivariant C S) (z : ι → Z) (R : ι → ℝ) {ε : ℝ}
    (hε : 0 ≤ ε) (k a : ℕ) :
    1 - ((k + a : ℝ) + 2 * (strayCount C (S z) R ε : ℝ) / (a + 1)) / C.card ≤
      (((permsOn C).filter fun σ : Equiv.Perm ι =>
        AcceptsPairwiseBy k C (S (z ∘ σ)) (2 * ε) t).card : ℝ) / (permsOn C).card := by
  have hc : (0 : ℝ) < C.card := by exact_mod_cast card_pos.2 ⟨t, ht⟩
  have hfilt : ((permsOn C).filter fun σ : Equiv.Perm ι =>
      AcceptsPairwiseBy k C (S (z ∘ σ)) (2 * ε) t) =
      (permsOn C).filter fun σ : Equiv.Perm ι => AcceptsPairwiseBy k C (S z) (2 * ε) (σ t) := by
    refine filter_congr fun σ hσ => ?_
    have : S (z ∘ σ) = fun i j => S z (σ i) (σ j) :=
      funext fun i => funext fun j => hS σ hσ z i j
    rw [this, acceptsPairwiseBy_comp hσ]
  rw [hfilt, ratio_permsOn ht (fun u => AcceptsPairwiseBy k C (S z) (2 * ε) u)]
  have hsplit := C.card_filter_add_card_filter_not
    (fun u => AcceptsPairwiseBy k C (S z) (2 * ε) u)
  have href := card_refused_by_le C (S z) R hε k a
  have e : ((C.filter fun u => AcceptsPairwiseBy k C (S z) (2 * ε) u).card : ℝ) =
      C.card - ((C.filter fun u => ¬ AcceptsPairwiseBy k C (S z) (2 * ε) u).card : ℝ) := by
    have : ((C.filter fun u => AcceptsPairwiseBy k C (S z) (2 * ε) u).card : ℝ) +
        ((C.filter fun u => ¬ AcceptsPairwiseBy k C (S z) (2 * ε) u).card : ℝ) = C.card := by
      exact_mod_cast hsplit
    linarith
  rw [e, sub_div, div_self hc.ne']
  have := div_le_div_of_nonneg_right href hc.le
  linarith

end Tournament

/-! ### Slice 14: shift -/

section Shift

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- Summing over the orderings what the test unit lands on: each unit of `C` is landed on by
`|orderings|/|C|` of them. -/
theorem card_mul_sum_permsOn {C : Finset ι} {t : ι} (ht : t ∈ C) (f : ι → ℝ) :
    (C.card : ℝ) * ∑ σ ∈ permsOn C, f (σ t) = (permsOn C).card * ∑ u ∈ C, f u := by
  have hfib : ∀ σ ∈ permsOn C, f (σ t) = ∑ u ∈ C, if σ t = u then f u else 0 := fun σ hσ => by
    rw [sum_ite_eq, if_pos ((mem_iff_of_mem_permsOn hσ t).2 ht)]
  rw [sum_congr rfl hfib, sum_comm, mul_sum, mul_sum]
  refine sum_congr rfl fun u hu => ?_
  rw [← sum_filter, sum_const, nsmul_eq_mul]
  have h := card_mul_card_filter_permsOn ht (fun v => v = u)
  simp only [filter_eq', if_pos hu, card_singleton, mul_one] at h
  have h' : (C.card : ℝ) * (((permsOn C).filter fun σ : Equiv.Perm ι => σ t = u).card : ℝ) =
      (permsOn C).card := by exact_mod_cast h
  calc (C.card : ℝ) * ((((permsOn C).filter fun σ : Equiv.Perm ι => σ t = u).card : ℝ) * f u)
      = ((C.card : ℝ) * (((permsOn C).filter fun σ : Equiv.Perm ι => σ t = u).card : ℝ)) * f u :=
        by ring
    _ = (permsOn C).card * f u := by rw [h']

omit [Fintype ι] in
/-- The weight of the units of `C` scoring at most `u`'s score, `u` included. The weighted rule
(Tibshirani, Barber, Candès and Ramdas 2019; Podkopaev and Ramdas 2021; Barber, Candès, Ramdas and
Tibshirani 2023) refuses `u` when this is at most `α`: their quantile with the candidate's mass at
`+∞`, in conformity form, ties counted for the candidate. -/
noncomputable def weightLe (C : Finset ι) (s q : ι → ℝ) (u : ι) : ℝ :=
  ∑ j ∈ C.filter (fun j => s j ≤ s u), q j

omit [Fintype ι] [DecidableEq ι] in
/-- **The weighted quantile lemma** (their Lemma 3; Barber et al.'s strange points): the units the
weighted rule refuses carry weight at most `α`, for any nonnegative weights. The refused units are
those scoring at most the highest refused one, whose weight-below is at most `α`. -/
theorem sum_weight_refused_le (C : Finset ι) (s q : ι → ℝ) (hq : ∀ j ∈ C, 0 ≤ q j) {α : ℝ}
    (hα : 0 ≤ α) : ∑ u ∈ C.filter (fun u => weightLe C s q u ≤ α), q u ≤ α := by
  rcases (C.filter fun u => weightLe C s q u ≤ α).eq_empty_or_nonempty with h | h
  · rw [h, sum_empty]
    exact hα
  · obtain ⟨u, hu, hmax⟩ := (C.filter fun u => weightLe C s q u ≤ α).exists_max_image s h
    have hsub : (C.filter fun u => weightLe C s q u ≤ α) ⊆ C.filter fun j => s j ≤ s u :=
      fun v hv => mem_filter.2 ⟨(mem_filter.1 hv).1, hmax v hv⟩
    calc ∑ v ∈ C.filter (fun u => weightLe C s q u ≤ α), q v
        ≤ ∑ j ∈ C.filter (fun j => s j ≤ s u), q j :=
          sum_le_sum_of_subset_of_nonneg hsub fun j hj _ => hq j (mem_filter.1 hj).1
      _ = weightLe C s q u := rfl
      _ ≤ α := (mem_filter.1 hu).2

omit [Fintype ι] [DecidableEq ι] in
theorem weightLe_const (C : Finset ι) (s : ι → ℝ) (c : ℝ) (u : ι) :
    weightLe C s (fun _ => c) u = rankLe C s u * c := by
  simp only [weightLe, sum_const, nsmul_eq_mul, rankLe]

omit [Fintype ι] [DecidableEq ι] in
/-- **At equal weights the weighted rule is CB-1**, at `k = ⌊α|C|⌋`. -/
theorem weighted_eq_cb1 (C : Finset ι) (s : ι → ℝ) {α : ℝ} (hα : 0 ≤ α) (hC : 0 < C.card)
    (u : ι) : α < weightLe C s (fun _ => (1 : ℝ) / C.card) u ↔ Accepts ⌊α * C.card⌋₊ C s u := by
  have hc : (0 : ℝ) < C.card := by exact_mod_cast hC
  rw [weightLe_const, Accepts, Nat.floor_lt (mul_nonneg hα (Nat.cast_nonneg _)), mul_one_div,
    lt_div_iff₀ hc]

theorem weightLe_comp {C : Finset ι} {σ : Equiv.Perm ι} (hσ : σ ∈ permsOn C) (s q : ι → ℝ)
    (u : ι) : weightLe C (s ∘ σ) (q ∘ σ) u = weightLe C s q (σ u) := by
  unfold weightLe
  exact sum_equiv σ (fun j => by simp only [mem_filter, Function.comp_apply,
    mem_iff_of_mem_permsOn hσ]) fun j _ => rfl

/-- **Weighted conformal under covariate or label shift**, per data set. The test unit's value
comes from the shifted law, with likelihood ratio `w` against the training law: `w(x)` for
covariate shift (Tibshirani et al., Corollary 1), `w(y)` for label shift (Podkopaev and Ramdas,
Theorem 2). Over the orderings of a data set the test unit lands on unit `u` with weight
`w(z_u)`. The weighted rule, with the units' normalized weights, refuses it with weight at most
`α` of the total. -/
theorem coverage_weighted {Z : Type*} {C : Finset ι} {t : ι} (ht : t ∈ C)
    {S : (ι → Z) → ι → ℝ} (hS : Equivariant C S) (z : ι → Z) (w : Z → ℝ) (hw : ∀ v, 0 ≤ w v)
    {α : ℝ} (hα : 0 ≤ α) :
    ∑ σ ∈ permsOn C, (if α < weightLe C (S (z ∘ σ))
        (fun j => w (z (σ j)) / ∑ i ∈ C, w (z i)) t then 0 else w (z (σ t))) ≤
      α * ∑ σ ∈ permsOn C, w (z (σ t)) := by
  have hc : (0 : ℝ) < C.card := by exact_mod_cast card_pos.2 ⟨t, ht⟩
  have hpt : ∀ σ ∈ permsOn C, (if α < weightLe C (S (z ∘ σ))
      (fun j => w (z (σ j)) / ∑ i ∈ C, w (z i)) t then 0 else w (z (σ t))) =
      (fun u => if α < weightLe C (S z) (fun j => w (z j) / ∑ i ∈ C, w (z i)) u then 0
        else w (z u)) (σ t) := by
    intro σ hσ
    have hfun : S (z ∘ σ) = S z ∘ σ := funext fun i => hS σ hσ z i
    rw [hfun, show (fun j => w (z (σ j)) / ∑ i ∈ C, w (z i)) =
      (fun j => w (z j) / ∑ i ∈ C, w (z i)) ∘ σ from rfl, weightLe_comp hσ]
  rw [sum_congr rfl hpt]
  have hkey : ∑ u ∈ C, (if α < weightLe C (S z) (fun j => w (z j) / ∑ i ∈ C, w (z i)) u then 0
      else w (z u)) ≤ α * ∑ i ∈ C, w (z i) := by
    have e : ∀ u ∈ C, (if α < weightLe C (S z) (fun j => w (z j) / ∑ i ∈ C, w (z i)) u then 0
        else w (z u)) = if weightLe C (S z) (fun j => w (z j) / ∑ i ∈ C, w (z i)) u ≤ α
          then w (z u) else 0 := fun u _ => by
      split_ifs with h1 h2 h2 <;> first | rfl | (exfalso; linarith)
    rw [sum_congr rfl e, ← sum_filter]
    rcases (sum_nonneg fun i (_ : i ∈ C) => hw (z i)).eq_or_lt with h0 | hpos
    · have hz : ∀ i ∈ C, w (z i) = 0 :=
        (sum_eq_zero_iff_of_nonneg fun i _ => hw (z i)).1 h0.symm
      rw [sum_eq_zero fun u hu => hz u (mem_filter.1 hu).1]
      exact le_of_eq (by rw [← h0, mul_zero])
    · have hne := hpos.ne'
      have hl := sum_weight_refused_le C (S z) (fun j => w (z j) / ∑ i ∈ C, w (z i))
        (fun j _ => div_nonneg (hw _) hpos.le) hα
      have e2 : ∑ u ∈ C.filter (fun u => weightLe C (S z)
          (fun j => w (z j) / ∑ i ∈ C, w (z i)) u ≤ α), w (z u) =
          (∑ i ∈ C, w (z i)) * ∑ u ∈ C.filter (fun u => weightLe C (S z)
            (fun j => w (z j) / ∑ i ∈ C, w (z i)) u ≤ α), w (z u) / ∑ i ∈ C, w (z i) := by
        rw [mul_sum]
        exact sum_congr rfl fun u _ => by field_simp
      rw [e2, mul_comm α]
      exact mul_le_mul_of_nonneg_left hl hpos.le
  refine (mul_le_mul_iff_of_pos_left hc).1 ?_
  rw [card_mul_sum_permsOn ht (fun u => if α < weightLe C (S z)
      (fun j => w (z j) / ∑ i ∈ C, w (z i)) u then 0 else w (z u)), mul_left_comm,
    card_mul_sum_permsOn ht (fun u => w (z u))]
  have hG : (0 : ℝ) ≤ (permsOn C).card := Nat.cast_nonneg _
  calc ((permsOn C).card : ℝ) * ∑ u ∈ C, (if α < weightLe C (S z)
        (fun j => w (z j) / ∑ i ∈ C, w (z i)) u then 0 else w (z u))
      ≤ (permsOn C).card * (α * ∑ i ∈ C, w (z i)) := mul_le_mul_of_nonneg_left hkey hG
    _ = α * ((permsOn C).card * ∑ u ∈ C, w (z u)) := by ring

omit [Fintype ι] [DecidableEq ι] in
/-- **The loss with estimated weights**, per data set: the units the rule refuses under estimated
weights `ph` carry true weight `p` at most `α` plus `Σ (p − ph)⁺`, the total variation between the
normalized true and estimated weights. Our statement: the papers give no finite bound with
estimated weights (Podkopaev and Ramdas's Corollary 1 is asymptotic). -/
theorem sum_weight_refused_estimated_le (C : Finset ι) (s p ph : ι → ℝ)
    (hph : ∀ j ∈ C, 0 ≤ ph j) {α : ℝ} (hα : 0 ≤ α) :
    ∑ u ∈ C.filter (fun u => weightLe C s ph u ≤ α), p u ≤
      α + ∑ u ∈ C, max 0 (p u - ph u) := by
  have h1 := sum_weight_refused_le C s ph hph hα
  have h2 : ∑ u ∈ C.filter (fun u => weightLe C s ph u ≤ α), p u ≤
      ∑ u ∈ C.filter (fun u => weightLe C s ph u ≤ α), ph u +
        ∑ u ∈ C.filter (fun u => weightLe C s ph u ≤ α), max 0 (p u - ph u) := by
    rw [← sum_add_distrib]
    exact sum_le_sum fun u _ => by linarith [le_max_right 0 (p u - ph u)]
  have h3 : ∑ u ∈ C.filter (fun u => weightLe C s ph u ≤ α), max 0 (p u - ph u) ≤
      ∑ u ∈ C, max 0 (p u - ph u) :=
    sum_le_sum_of_subset_of_nonneg (filter_subset _ _) fun _ _ _ => le_max_left _ _
  linarith

/-- Swapping a unit with the test position moves weight between them only. -/
theorem weightLe_swap_le (s q : ι → ℝ) (K t : ι) (hq : q K ≤ q t) :
    weightLe univ (s ∘ Equiv.swap K t) q K ≤ weightLe univ s q t := by
  have hcomp : weightLe univ (s ∘ Equiv.swap K t) q K =
      weightLe univ s (q ∘ Equiv.swap K t) t := by
    have h := weightLe_comp (mem_permsOn_univ (Equiv.swap K t)) s (q ∘ Equiv.swap K t) K
    simp only [Equiv.swap_apply_left] at h
    rw [← h]
    congr 1
    funext j
    simp [Function.comp_apply, Equiv.swap_apply_self]
  rw [hcomp]
  by_cases hKt : K = t
  · subst hKt
    simp [Equiv.swap_self]
  unfold weightLe
  rw [sum_filter, sum_filter, ← sub_nonneg, ← sum_sub_distrib,
    Fintype.sum_eq_add K t hKt fun j hj => by
      simp [Equiv.swap_apply_of_ne_of_ne hj.1 hj.2]]
  simp only [Function.comp_apply, Equiv.swap_apply_left, Equiv.swap_apply_right, le_refl,
    ↓reduceIte]
  split_ifs <;> linarith

/-- **Beyond exchangeability** (Barber, Candès, Ramdas and Tibshirani 2023, Theorem 2), on a
finite law. `P` is any law of the data, the score symmetric in the units (`Equivariant univ`), and
the weights fixed by position, nonnegative and summing to one, the test position's the largest.
Then the weighted rule refuses the test unit with chance at most `α` plus
`Σ_K q_K · TV(P, P ∘ swap(K, t))`. With `P` exchangeable the swaps cost nothing: their
coverage gap is the weighted swap distances. -/
theorem refused_le_swap_tv {Z : Type*} [Fintype Z] {S : (ι → Z) → ι → ℝ}
    (hS : Equivariant univ S) (P : (ι → Z) → ℝ) (hP0 : ∀ z, 0 ≤ P z) (q : ι → ℝ)
    (hq0 : ∀ i, 0 ≤ q i) (hq1 : ∑ i, q i = 1) (t : ι) (hqt : ∀ i, q i ≤ q t) {α : ℝ}
    (hα : 0 ≤ α) :
    ∑ z, P z * (if weightLe univ (S z) q t ≤ α then 1 else 0) ≤
      α * ∑ z, P z + ∑ K, q K * ∑ z, max 0 (P (z ∘ Equiv.swap K t) - P z) := by
  have hstep : ∀ z K, weightLe univ (S z) q t ≤ α →
      weightLe univ (S (z ∘ Equiv.swap K t)) q K ≤ α := by
    intro z K h
    have hfun : S (z ∘ Equiv.swap K t) = S z ∘ Equiv.swap K t :=
      funext fun i => hS _ (mem_permsOn_univ _) z i
    rw [hfun]
    exact (weightLe_swap_le (S z) q K t (hqt K)).trans h
  have hinv : ∀ K (z : ι → Z), (z ∘ Equiv.swap K t) ∘ Equiv.swap K t = z := fun K z => by
    funext j
    simp [Function.comp_apply, Equiv.swap_apply_self]
  have hK : ∀ K, ∑ z, P z * (if weightLe univ (S (z ∘ Equiv.swap K t)) q K ≤ α then 1 else 0) ≤
      ∑ z, P z * (if weightLe univ (S z) q K ≤ α then 1 else 0) +
        ∑ z, max 0 (P (z ∘ Equiv.swap K t) - P z) := by
    intro K
    have hre := sum_comp_perm (fun z => P (z ∘ Equiv.swap K t) *
      (if weightLe univ (S z) q K ≤ α then (1 : ℝ) else 0)) (Equiv.swap K t)
    simp only [hinv] at hre
    rw [hre, ← sum_add_distrib]
    refine sum_le_sum fun z _ => ?_
    split_ifs
    · linarith [le_max_right 0 (P (z ∘ Equiv.swap K t) - P z)]
    · linarith [le_max_left 0 (P (z ∘ Equiv.swap K t) - P z)]
  have hlem : ∀ z, ∑ K, q K * (if weightLe univ (S z) q K ≤ α then 1 else 0) ≤ α := by
    intro z
    have := sum_weight_refused_le univ (S z) q (fun j _ => hq0 j) hα
    rw [sum_filter] at this
    simpa only [mul_ite, mul_one, mul_zero] using this
  calc ∑ z, P z * (if weightLe univ (S z) q t ≤ α then 1 else 0)
      = ∑ K, q K * ∑ z, P z * (if weightLe univ (S z) q t ≤ α then 1 else 0) := by
        rw [← sum_mul, hq1, one_mul]
    _ ≤ ∑ K, q K * ∑ z, P z *
          (if weightLe univ (S (z ∘ Equiv.swap K t)) q K ≤ α then 1 else 0) := by
        refine sum_le_sum fun K _ => mul_le_mul_of_nonneg_left
          (sum_le_sum fun z _ => mul_le_mul_of_nonneg_left ?_ (hP0 z)) (hq0 K)
        split_ifs with h1 h2 <;> first | exact absurd (hstep z K h1) h2 | norm_num
    _ ≤ ∑ K, q K * (∑ z, P z * (if weightLe univ (S z) q K ≤ α then 1 else 0) +
          ∑ z, max 0 (P (z ∘ Equiv.swap K t) - P z)) :=
        sum_le_sum fun K _ => mul_le_mul_of_nonneg_left (hK K) (hq0 K)
    _ = ∑ z, P z * ∑ K, q K * (if weightLe univ (S z) q K ≤ α then 1 else 0) +
          ∑ K, q K * ∑ z, max 0 (P (z ∘ Equiv.swap K t) - P z) := by
        simp only [mul_add, sum_add_distrib]
        congr 1
        simp only [mul_sum]
        rw [sum_comm]
        exact sum_congr rfl fun z _ => sum_congr rfl fun K _ => by ring
    _ ≤ ∑ z, P z * α + ∑ K, q K * ∑ z, max 0 (P (z ∘ Equiv.swap K t) - P z) := by
        exact add_le_add (sum_le_sum fun z _ => mul_le_mul_of_nonneg_left (hlem z) (hP0 z)) le_rfl
    _ = α * ∑ z, P z + ∑ K, q K * ∑ z, max 0 (P (z ∘ Equiv.swap K t) - P z) := by
        rw [mul_sum]
        congr 1
        exact sum_congr rfl fun z _ => by ring

end Shift

/-! ### Slice 14's follow-up 1: the weighted predictive distribution's PIT

The weighted predictive distribution read at the test unit's value, with the tied mass randomized,
is uniform given the data set; read with ties counted for the candidate it is `weightLe`, whose
super-uniformity at every level is `coverage_weighted`. -/

section WeightedPit

/-- Clamps add along a split of `[0, Q + b]` at `Q`. -/
theorem clamp_add {Q b x : ℝ} (hQ : 0 ≤ Q) (hb : 0 ≤ b) :
    max 0 (min Q x) + max 0 (min b (x - Q)) = max 0 (min (Q + b) x) := by
  rcases le_total x 0 with h1 | h1
  · rw [min_eq_right (h1.trans hQ), max_eq_left h1, min_eq_right (by linarith : x - Q ≤ b),
      max_eq_left (by linarith : x - Q ≤ 0), min_eq_right (by linarith : x ≤ Q + b),
      max_eq_left h1, add_zero]
  rcases le_total x Q with h2 | h2
  · rw [min_eq_right h2, max_eq_right h1, min_eq_right (by linarith : x - Q ≤ b),
      max_eq_left (by linarith : x - Q ≤ 0), min_eq_right (by linarith : x ≤ Q + b),
      max_eq_right h1, add_zero]
  rcases le_total x (Q + b) with h3 | h3
  · rw [min_eq_left h2, max_eq_right hQ, min_eq_right (by linarith : x - Q ≤ b),
      max_eq_right (by linarith : 0 ≤ x - Q), min_eq_right h3, max_eq_right h1]
    ring
  · rw [min_eq_left h2, max_eq_right hQ, min_eq_left (by linarith : b ≤ x - Q), max_eq_right hb,
      min_eq_left h3, max_eq_right (by linarith : 0 ≤ Q + b)]

/-- **Levels tile** `[0, Σ b]`: level `v` covers `[a_v, a_v + b_v)`, `a_v` the mass of the levels
below it, so the clamps of `x` into the levels sum to `x` clamped into `[0, Σ b]`. -/
theorem levels_tile (T : Finset ℝ) (b : ℝ → ℝ) (hb : ∀ v ∈ T, 0 ≤ b v) (x : ℝ) :
    ∑ v ∈ T, max 0 (min (b v) (x - ∑ w ∈ T.filter (· < v), b w)) =
      max 0 (min (∑ v ∈ T, b v) x) := by
  induction T using Finset.induction_on_max with
  | empty => simp
  | insert M T hM ih =>
    have hMT : M ∉ T := fun h => lt_irrefl M (hM M h)
    have hbT : ∀ v ∈ T, 0 ≤ b v := fun v hv => hb v (mem_insert_of_mem hv)
    have hM' : (insert M T).filter (· < M) = T := by
      rw [filter_insert, if_neg (lt_irrefl M), filter_true_of_mem hM]
    have hrest : ∑ v ∈ T, max 0 (min (b v) (x - ∑ w ∈ (insert M T).filter (· < v), b w)) =
        ∑ v ∈ T, max 0 (min (b v) (x - ∑ w ∈ T.filter (· < v), b w)) :=
      sum_congr rfl fun v hv => by rw [filter_insert, if_neg (not_lt.2 (hM v hv).le)]
    rw [sum_insert hMT, hM', hrest, ih hbT, sum_insert hMT, add_comm (b M),
      add_comm (max 0 (min (b M) _)), clamp_add (sum_nonneg hbT) (hb M (mem_insert_self M T))]

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

omit [Fintype ι] in
/-- The weight of the units of `C` scoring strictly below `u`'s score. -/
noncomputable def weightLt (C : Finset ι) (s q : ι → ℝ) (u : ι) : ℝ :=
  ∑ j ∈ C.filter (fun j => s j < s u), q j

omit [Fintype ι] in
/-- The weight of the units of `C` tied with `u`, `u` included. -/
noncomputable def weightEq (C : Finset ι) (s q : ι → ℝ) (u : ι) : ℝ :=
  ∑ j ∈ C.filter (fun j => s j = s u), q j

omit [Fintype ι] in
/-- **The weighted predictive distribution's PIT, its CDF**: `P(T ≤ x)` for `T = A + V B`, with
`A` the weight scoring strictly below `u` (`weightLt`), `B` the weight tied with it (`weightEq`),
and `V` uniform on `(0, 1]`. At `V = 1` the PIT is `weightLe`. -/
noncomputable def wPitCdf (C : Finset ι) (s q : ι → ℝ) (u : ι) (x : ℝ) : ℝ :=
  if weightEq C s q u = 0 then (if weightLt C s q u ≤ x then 1 else 0)
  else unitClamp ((x - weightLt C s q u) / weightEq C s q u)

omit [Fintype ι] [DecidableEq ι] in
theorem wPitCdf_nonneg (C : Finset ι) (s q : ι → ℝ) (u : ι) (x : ℝ) : 0 ≤ wPitCdf C s q u x := by
  unfold wPitCdf unitClamp
  split_ifs <;> norm_num

omit [Fintype ι] [DecidableEq ι] in
theorem wPitCdf_le_one (C : Finset ι) (s q : ι → ℝ) (u : ι) (x : ℝ) : wPitCdf C s q u x ≤ 1 := by
  unfold wPitCdf unitClamp
  split_ifs <;> norm_num

theorem weightLt_comp {C : Finset ι} {σ : Equiv.Perm ι} (hσ : σ ∈ permsOn C) (s q : ι → ℝ)
    (u : ι) : weightLt C (s ∘ σ) (q ∘ σ) u = weightLt C s q (σ u) := by
  unfold weightLt
  exact sum_equiv σ (fun j => by simp only [mem_filter, Function.comp_apply,
    mem_iff_of_mem_permsOn hσ]) fun j _ => rfl

theorem weightEq_comp {C : Finset ι} {σ : Equiv.Perm ι} (hσ : σ ∈ permsOn C) (s q : ι → ℝ)
    (u : ι) : weightEq C (s ∘ σ) (q ∘ σ) u = weightEq C s q (σ u) := by
  unfold weightEq
  exact sum_equiv σ (fun j => by simp only [mem_filter, Function.comp_apply,
    mem_iff_of_mem_permsOn hσ]) fun j _ => rfl

theorem wPitCdf_comp {C : Finset ι} {σ : Equiv.Perm ι} (hσ : σ ∈ permsOn C) (s q : ι → ℝ)
    (u : ι) (x : ℝ) : wPitCdf C (s ∘ σ) (q ∘ σ) u x = wPitCdf C s q (σ u) x := by
  unfold wPitCdf
  rw [weightLt_comp hσ, weightEq_comp hσ]

/-- `B · P(V ≤ y/B) = y` clamped into `[0, B]`. -/
theorem mul_unitClamp_div {B y : ℝ} (hB : 0 < B) : B * unitClamp (y / B) = max 0 (min B y) := by
  unfold unitClamp
  have h : B * (y / B) = y := by field_simp
  rw [mul_max_of_nonneg _ _ hB.le, mul_min_of_nonneg _ _ hB.le, mul_zero, mul_one, h]

omit [Fintype ι] [DecidableEq ι] in
/-- **The weighted PIT is uniform given the data set**: with the test unit on unit `u` with chance
`q u`, `P(T ≤ x) = Σ_u q_u P(A_u + V B_u ≤ x)` is `x` clamped into `[0, Σ q]`. The tie levels'
intervals `[A, A + B)` tile `[0, Σ q]`. -/
theorem sum_mul_wPitCdf (C : Finset ι) (s q : ι → ℝ) (hq : ∀ j ∈ C, 0 ≤ q j) (x : ℝ) :
    ∑ u ∈ C, q u * wPitCdf C s q u x = max 0 (min (∑ u ∈ C, q u) x) := by
  set Bv : ℝ → ℝ := fun v => ∑ j ∈ C.filter (fun j => s j = v), q j with hBv
  set Av : ℝ → ℝ := fun v => ∑ j ∈ C.filter (fun j => s j < v), q j with hAv
  set g : ℝ → ℝ := fun v => max 0 (min (Bv v) (x - Av v)) with hg
  have hBnn : ∀ v, 0 ≤ Bv v := fun v => sum_nonneg fun j hj => hq j (mem_filter.1 hj).1
  -- each unit: its chance times its weight is its share of its level's clamp
  have hunit : ∀ u ∈ C, q u * wPitCdf C s q u x = q u / Bv (s u) * g (s u) := by
    intro u hu
    have hle : q u ≤ Bv (s u) :=
      single_le_sum (fun j hj => hq j (mem_filter.1 hj).1) (mem_filter.2 ⟨hu, rfl⟩)
    unfold wPitCdf
    change q u * (if Bv (s u) = 0 then (if Av (s u) ≤ x then 1 else 0)
      else unitClamp ((x - Av (s u)) / Bv (s u))) = q u / Bv (s u) * g (s u)
    by_cases h0 : Bv (s u) = 0
    · have : q u = 0 := le_antisymm (h0 ▸ hle) (hq u hu)
      rw [this, zero_mul, zero_div, zero_mul]
    · have hpos : 0 < Bv (s u) := lt_of_le_of_ne (hBnn _) (Ne.symm h0)
      rw [if_neg h0]
      change q u * unitClamp ((x - Av (s u)) / Bv (s u)) =
        q u / Bv (s u) * max 0 (min (Bv (s u)) (x - Av (s u)))
      rw [← mul_unitClamp_div hpos, ← mul_assoc, div_mul_cancel₀ _ h0]
  rw [sum_congr rfl hunit]
  -- units to levels
  have hlev : ∑ u ∈ C, q u / Bv (s u) * g (s u) = ∑ v ∈ C.image s, g v := by
    refine (sum_image' (fun u => q u / Bv (s u) * g (s u)) fun i hi => ?_).symm
    have e : ∑ j ∈ C with s j = s i, q j / Bv (s j) * g (s j) =
        (∑ j ∈ C with s j = s i, q j) / Bv (s i) * g (s i) := by
      rw [sum_div, sum_mul]
      exact sum_congr rfl fun j hj => by rw [(mem_filter.1 hj).2]
    rw [e]
    change g (s i) = Bv (s i) / Bv (s i) * g (s i)
    rcases eq_or_ne (Bv (s i)) 0 with h0 | h0
    · rw [h0, div_zero, zero_mul]
      change max 0 (min (Bv (s i)) (x - Av (s i))) = 0
      rw [h0]
      exact max_eq_left (min_le_left _ _)
    · rw [div_self h0, one_mul]
  rw [hlev]
  -- the mass below a level is the levels below it
  have hA : ∀ v, Av v = ∑ w ∈ (C.image s).filter (· < v), Bv w := by
    intro v
    rw [show (C.image s).filter (· < v) = (C.filter fun j => s j < v).image s by
      ext w
      simp only [mem_filter, mem_image]
      constructor
      · rintro ⟨⟨j, hj, rfl⟩, hw⟩
        exact ⟨j, ⟨hj, hw⟩, rfl⟩
      · rintro ⟨j, ⟨hj, hw⟩, rfl⟩
        exact ⟨⟨j, hj, rfl⟩, hw⟩]
    refine (sum_image' q fun i hi => ?_).symm
    have hiv : s i < v := (mem_filter.1 hi).2
    rw [filter_filter]
    exact sum_congr (filter_congr fun j _ => ⟨fun h => ⟨by rw [h]; exact hiv, h⟩,
      fun h => h.2⟩) fun _ _ => rfl
  have hsumB : ∑ v ∈ C.image s, Bv v = ∑ u ∈ C, q u :=
    sum_image' q fun i _ => rfl
  rw [← hsumB, ← levels_tile (C.image s) Bv (fun v _ => hBnn v) x]
  refine sum_congr rfl fun v _ => ?_
  change max 0 (min (Bv v) (x - Av v)) = _
  rw [hA v]

/-- **The weighted predictive distribution's PIT, over orderings** (Tibshirani, Barber, Candès and
Ramdas 2019 for covariate shift, `w(x)`; Podkopaev and Ramdas 2021 for response shift, `w(y)`).
Over the orderings of a data set the test unit lands on unit `u` with weight `w(z_u)`; read at its
value with the tied mass randomized, the weighted distribution's PIT is exactly uniform. The score
is equivariant, so this is the full form; the split form is `S z = s ∘ z`. -/
theorem pit_weighted {Z : Type*} {C : Finset ι} {t : ι} (ht : t ∈ C)
    {S : (ι → Z) → ι → ℝ} (hS : Equivariant C S) (z : ι → Z) (w : Z → ℝ) (hw : ∀ v, 0 ≤ w v)
    {x : ℝ} (hx0 : 0 ≤ x) (hx1 : x ≤ 1) :
    ∑ σ ∈ permsOn C, w (z (σ t)) *
        wPitCdf C (S (z ∘ σ)) (fun j => w (z (σ j)) / ∑ i ∈ C, w (z i)) t x =
      x * ∑ σ ∈ permsOn C, w (z (σ t)) := by
  set q : ι → ℝ := fun j => w (z j) / ∑ i ∈ C, w (z i) with hq
  have hpt : ∀ σ ∈ permsOn C, w (z (σ t)) *
      wPitCdf C (S (z ∘ σ)) (fun j => w (z (σ j)) / ∑ i ∈ C, w (z i)) t x =
      (fun u => w (z u) * wPitCdf C (S z) q u x) (σ t) := by
    intro σ hσ
    have hfun : S (z ∘ σ) = S z ∘ σ := funext fun i => hS σ hσ z i
    rw [hfun, show (fun j => w (z (σ j)) / ∑ i ∈ C, w (z i)) = q ∘ σ from rfl, wPitCdf_comp hσ]
  rw [sum_congr rfl hpt]
  have hc : (0 : ℝ) < C.card := by exact_mod_cast card_pos.2 ⟨t, ht⟩
  have hkey : ∑ u ∈ C, w (z u) * wPitCdf C (S z) q u x = x * ∑ i ∈ C, w (z i) := by
    rcases (sum_nonneg fun i (_ : i ∈ C) => hw (z i)).eq_or_lt with h0 | hpos
    · have hz : ∀ i ∈ C, w (z i) = 0 :=
        (sum_eq_zero_iff_of_nonneg fun i _ => hw (z i)).1 h0.symm
      rw [sum_eq_zero fun u hu => by rw [hz u hu, zero_mul], ← h0, mul_zero]
    · have hW := hpos.ne'
      have htile := sum_mul_wPitCdf C (S z) q (fun j _ => div_nonneg (hw _) hpos.le) x
      have hq1 : ∑ u ∈ C, q u = 1 := by
        change ∑ u ∈ C, w (z u) / ∑ i ∈ C, w (z i) = 1
        rw [← sum_div, div_self hW]
      rw [hq1, min_eq_right hx1, max_eq_right hx0] at htile
      calc ∑ u ∈ C, w (z u) * wPitCdf C (S z) q u x =
            (∑ i ∈ C, w (z i)) * ∑ u ∈ C, q u * wPitCdf C (S z) q u x := by
            rw [mul_sum]
            refine sum_congr rfl fun u _ => ?_
            change w (z u) * wPitCdf C (S z) q u x =
              (∑ i ∈ C, w (z i)) * ((w (z u) / ∑ i ∈ C, w (z i)) * wPitCdf C (S z) q u x)
            rw [← mul_assoc, mul_comm (∑ i ∈ C, w (z i)), div_mul_cancel₀ _ hW]
        _ = x * ∑ i ∈ C, w (z i) := by rw [htile, mul_comm]
  refine mul_left_cancel₀ hc.ne' ?_
  rw [card_mul_sum_permsOn ht (fun u => w (z u) * wPitCdf C (S z) q u x), hkey,
    mul_left_comm (C.card : ℝ) x, card_mul_sum_permsOn ht (fun u => w (z u))]
  ring

omit [Fintype ι] [DecidableEq ι] in
/-- **The PIT with an estimated ratio** (our statement): with the test unit landing by the true
normalized weights `p` and the PIT read with estimates `ph`, `P(T ≤ x)` is within the total
variation `Σ (p − ph)⁺` of `x`, from either side. -/
theorem wPit_estimated (C : Finset ι) (s p ph : ι → ℝ) (hph : ∀ j ∈ C, 0 ≤ ph j)
    (hp1 : ∑ j ∈ C, p j = 1) (hph1 : ∑ j ∈ C, ph j = 1) {x : ℝ} (hx0 : 0 ≤ x) (hx1 : x ≤ 1) :
    |∑ u ∈ C, p u * wPitCdf C s ph u x - x| ≤ ∑ u ∈ C, max 0 (p u - ph u) := by
  have htile := sum_mul_wPitCdf C s ph hph x
  rw [hph1, min_eq_right hx1, max_eq_right hx0] at htile
  have hd : ∑ u ∈ C, p u * wPitCdf C s ph u x - x =
      ∑ u ∈ C, (p u - ph u) * wPitCdf C s ph u x := by
    have e : ∑ u ∈ C, p u * wPitCdf C s ph u x - x =
        ∑ u ∈ C, p u * wPitCdf C s ph u x - ∑ u ∈ C, ph u * wPitCdf C s ph u x := by
      rw [htile]
    rw [e, ← sum_sub_distrib]
    exact sum_congr rfl fun u _ => by ring
  have hterm : ∀ u, max 0 (p u - ph u) - max 0 (ph u - p u) = p u - ph u := fun u => by
    rcases le_total 0 (p u - ph u) with h | h
    · rw [max_eq_right h, max_eq_left (by linarith : ph u - p u ≤ 0)]
      ring
    · rw [max_eq_left h, max_eq_right (by linarith : 0 ≤ ph u - p u)]
      ring
  have hsym : ∑ u ∈ C, max 0 (ph u - p u) = ∑ u ∈ C, max 0 (p u - ph u) := by
    have h : ∑ u ∈ C, (max 0 (p u - ph u) - max 0 (ph u - p u)) = ∑ u ∈ C, (p u - ph u) :=
      sum_congr rfl fun u _ => hterm u
    rw [sum_sub_distrib, sum_sub_distrib, hp1, hph1] at h
    linarith
  rw [hd, abs_le]
  constructor
  · rw [← hsym, ← sum_neg_distrib]
    refine sum_le_sum fun u _ => ?_
    have h0 := wPitCdf_nonneg C s ph u x
    have h1 := wPitCdf_le_one C s ph u x
    rcases le_total 0 (p u - ph u) with h | h
    · rw [max_eq_left (by linarith : ph u - p u ≤ 0)]
      nlinarith
    · rw [max_eq_right (by linarith : 0 ≤ ph u - p u)]
      nlinarith
  · refine sum_le_sum fun u _ => ?_
    have h0 := wPitCdf_nonneg C s ph u x
    have h1 := wPitCdf_le_one C s ph u x
    rcases le_total 0 (p u - ph u) with h | h
    · rw [max_eq_right h]
      nlinarith
    · rw [max_eq_left h]
      nlinarith

end WeightedPit

/-! ### Slice 14's follow-up 2: the cross-fitted forms under weights

Barber, Candès, Ramdas and Tibshirani 2023: the nonexchangeable jackknife+ (their Theorem 5),
its `n`-fold cross-conformal form (their Appendix D), and Lei's lemma on weighted strange points
(their Lemma 2). The tournament correction is the cross-fitted rule with its own pair score, so it
follows. -/

section WeightedCrossFit

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

omit [Fintype ι] [DecidableEq ι] in
/-- The weight of the units of `C` that `u` is stranger than. The weighted cross-fitted rule refuses
`u` when this is at least `1 − α` of the total: their weighted strange points. -/
def strangeWeight {β : Type*} [LinearOrder β] (C : Finset ι) (s : ι → ι → β) (q : ι → ℝ)
    (u : ι) : ℝ :=
  ∑ j ∈ strangerThan C s u, q j

omit [Fintype ι] [DecidableEq ι] in
/-- **Weighted strange points** (Lei's lemma, Barber et al.'s Lemma 2): whatever the comparisons,
the units stranger than at least `1 − α` of the weight carry at most `2α` of it. Each refused unit
is stranger than all but `α` of the weight; among the refused units each comparison goes at most one
way, so they are stranger than at most half of their own weight's square. -/
theorem sum_weight_strange_le {β : Type*} [LinearOrder β] (C : Finset ι) (s : ι → ι → β)
    (q : ι → ℝ) (hq : ∀ j ∈ C, 0 ≤ q j) {α : ℝ} (hα : 0 ≤ α) :
    ∑ u ∈ C.filter (fun u => (1 - α) * ∑ j ∈ C, q j ≤ strangeWeight C s q u), q u ≤
      2 * α * ∑ j ∈ C, q j := by
  classical
  set F := C.filter (fun u => (1 - α) * ∑ j ∈ C, q j ≤ strangeWeight C s q u) with hF
  set W := ∑ j ∈ C, q j with hW
  set QF := ∑ u ∈ F, q u with hQF
  have hFC : F ⊆ C := filter_subset _ _
  have hqF : ∀ j ∈ F, 0 ≤ q j := fun j hj => hq j (hFC hj)
  set a : ι → ℝ := fun u => ∑ j ∈ F, if s u j < s j u then q j else 0 with ha
  have hunit : ∀ u ∈ F, (1 - α) * W ≤ a u + (W - QF) := by
    intro u hu
    have h1 := (mem_filter.1 hu).2
    have h2 : strangeWeight C s q u ≤ a u + (W - QF) := by
      have hsd := sum_sdiff hFC (f := fun j => if s u j < s j u then q j else 0)
      have hrest : ∑ j ∈ C \ F, (if s u j < s j u then q j else 0) ≤ ∑ j ∈ C \ F, q j :=
        sum_le_sum fun j hj => by
          split_ifs
          · exact le_rfl
          · exact hq j (mem_sdiff.1 hj).1
      have hW' : ∑ j ∈ C \ F, q j = W - QF := by
        rw [sum_sdiff_eq_sub hFC]
      unfold strangeWeight strangerThan
      rw [sum_filter, ← hsd]
      change ∑ j ∈ C \ F, (if s u j < s j u then q j else 0) + a u ≤ a u + (W - QF)
      linarith
    linarith
  have hsum : (1 - α) * W * QF ≤ ∑ u ∈ F, q u * a u + (W - QF) * QF := by
    have h := sum_le_sum fun u hu => mul_le_mul_of_nonneg_left (hunit u hu) (hqF u hu)
    rw [← sum_mul, mul_comm, ← hQF] at h
    simp only [mul_add, sum_add_distrib, ← sum_mul, ← hQF] at h
    linarith
  have hpair : 2 * ∑ u ∈ F, q u * a u ≤ QF * QF := by
    have hX : ∑ u ∈ F, q u * a u =
        ∑ u ∈ F, ∑ j ∈ F, if s u j < s j u then q u * q j else 0 := by
      refine sum_congr rfl fun u _ => ?_
      rw [ha, mul_sum]
      exact sum_congr rfl fun j _ => by split_ifs <;> ring
    have hY : ∑ u ∈ F, q u * a u =
        ∑ u ∈ F, ∑ j ∈ F, if s j u < s u j then q u * q j else 0 := by
      rw [hX, sum_comm]
      exact sum_congr rfl fun u _ => sum_congr rfl fun j _ => by split_ifs <;> ring
    have hsq : QF * QF = ∑ u ∈ F, ∑ j ∈ F, q u * q j := by
      rw [hQF, sum_mul_sum]
    rw [hsq, two_mul]
    nth_rw 1 [hX]
    rw [hY, ← sum_add_distrib]
    refine sum_le_sum fun u hu => ?_
    rw [← sum_add_distrib]
    refine sum_le_sum fun j hj => ?_
    have hprod : 0 ≤ q u * q j := mul_nonneg (hqF u hu) (hqF j hj)
    by_cases h1 : s u j < s j u
    · have h2 : ¬ s j u < s u j := fun h2 => lt_asymm h1 h2
      rw [if_pos h1, if_neg h2, add_zero]
    · rw [if_neg h1, zero_add]
      split_ifs
      · exact le_rfl
      · exact hprod
  have hQF0 : 0 ≤ QF := sum_nonneg hqF
  have hkey : QF * QF ≤ 2 * α * W * QF := by nlinarith
  rcases hQF0.eq_or_lt with h0 | hpos
  · rw [← h0]
    have hW0 : 0 ≤ W := sum_nonneg hq
    nlinarith
  · exact le_of_mul_le_mul_right (by linarith) hpos

/-- **The swap step**: with the test position's weight the largest, a unit `t` stranger than a
share of the weight, read after swapping the data at `K` and `t`, is a unit `K` stranger than at
least that share. -/
theorem strangeWeight_swap_ge {β : Type*} [LinearOrder β] (s : ι → ι → β) (q : ι → ℝ)
    (K t : ι) (hqt : q K ≤ q t) :
    strangeWeight univ s q t ≤
      strangeWeight univ (fun i j => s (Equiv.swap K t i) (Equiv.swap K t j)) q K := by
  unfold strangeWeight strangerThan
  rw [sum_filter, sum_filter]
  rw [← Equiv.sum_comp (Equiv.swap K t)
    (fun j => if s (Equiv.swap K t K) (Equiv.swap K t j) < s (Equiv.swap K t j)
      (Equiv.swap K t K) then q j else 0)]
  simp only [Equiv.swap_apply_left, Equiv.swap_apply_self]
  refine sum_le_sum fun j _ => ?_
  split_ifs with h
  · by_cases hjt : j = t
    · subst hjt
      exact absurd h (lt_irrefl _)
    by_cases hjK : j = K
    · subst hjK
      rw [Equiv.swap_apply_left]
      exact hqt
    · rw [Equiv.swap_apply_of_ne_of_ne hjK hjt]
  · exact le_rfl

/-- **The cross-fitted rule under weights, beyond exchangeability** (Barber, Candès, Ramdas and
Tibshirani 2023, Theorem 5 and Appendix D), on a finite law. `P` is any law of the data, the pair
score equivariant (a symmetric fit's leave-two-out score), the weights fixed by position, summing to
one, the test position's the largest. The rule refuses the test unit with chance at most `2α` plus
`Σ_K q_K · TV(P, P ∘ swap(K, t))`. -/
theorem crossFit_weighted_swap_tv {Z β : Type*} [Fintype Z] [LinearOrder β]
    {S : (ι → Z) → ι → ι → β} (hS : PairEquivariant univ S) (P : (ι → Z) → ℝ)
    (hP0 : ∀ z, 0 ≤ P z) (q : ι → ℝ) (hq0 : ∀ i, 0 ≤ q i) (hq1 : ∑ i, q i = 1) (t : ι)
    (hqt : ∀ i, q i ≤ q t) {α : ℝ} (hα : 0 ≤ α) :
    ∑ z, P z * (if 1 - α ≤ strangeWeight univ (S z) q t then 1 else 0) ≤
      2 * α * ∑ z, P z + ∑ K, q K * ∑ z, max 0 (P (z ∘ Equiv.swap K t) - P z) := by
  have hstep : ∀ z K, 1 - α ≤ strangeWeight univ (S z) q t →
      1 - α ≤ strangeWeight univ (S (z ∘ Equiv.swap K t)) q K := by
    intro z K h
    have hfun : S (z ∘ Equiv.swap K t) =
        fun i j => S z (Equiv.swap K t i) (Equiv.swap K t j) :=
      funext fun i => funext fun j => hS _ (mem_permsOn_univ _) z i j
    rw [hfun]
    exact h.trans (strangeWeight_swap_ge (S z) q K t (hqt K))
  have hinv : ∀ K (z : ι → Z), (z ∘ Equiv.swap K t) ∘ Equiv.swap K t = z := fun K z => by
    funext j
    simp [Function.comp_apply, Equiv.swap_apply_self]
  have hK : ∀ K, ∑ z, P z * (if 1 - α ≤ strangeWeight univ (S (z ∘ Equiv.swap K t)) q K
      then 1 else 0) ≤
      ∑ z, P z * (if 1 - α ≤ strangeWeight univ (S z) q K then 1 else 0) +
        ∑ z, max 0 (P (z ∘ Equiv.swap K t) - P z) := by
    intro K
    have hre := sum_comp_perm (fun z => P (z ∘ Equiv.swap K t) *
      (if 1 - α ≤ strangeWeight univ (S z) q K then (1 : ℝ) else 0)) (Equiv.swap K t)
    simp only [hinv] at hre
    rw [hre, ← sum_add_distrib]
    refine sum_le_sum fun z _ => ?_
    split_ifs
    · linarith [le_max_right 0 (P (z ∘ Equiv.swap K t) - P z)]
    · linarith [le_max_left 0 (P (z ∘ Equiv.swap K t) - P z)]
  have hlem : ∀ z, ∑ K, q K * (if 1 - α ≤ strangeWeight univ (S z) q K then 1 else 0) ≤
      2 * α := by
    intro z
    have := sum_weight_strange_le univ (S z) q (fun j _ => hq0 j) hα
    rw [hq1, mul_one, sum_filter] at this
    simpa only [mul_ite, mul_one, mul_zero] using this
  calc ∑ z, P z * (if 1 - α ≤ strangeWeight univ (S z) q t then 1 else 0)
      = ∑ K, q K * ∑ z, P z * (if 1 - α ≤ strangeWeight univ (S z) q t then 1 else 0) := by
        rw [← sum_mul, hq1, one_mul]
    _ ≤ ∑ K, q K * ∑ z, P z *
          (if 1 - α ≤ strangeWeight univ (S (z ∘ Equiv.swap K t)) q K then 1 else 0) := by
        refine sum_le_sum fun K _ => mul_le_mul_of_nonneg_left
          (sum_le_sum fun z _ => mul_le_mul_of_nonneg_left ?_ (hP0 z)) (hq0 K)
        split_ifs with h1 h2 <;> first | exact absurd (hstep z K h1) h2 | norm_num
    _ ≤ ∑ K, q K * (∑ z, P z * (if 1 - α ≤ strangeWeight univ (S z) q K then 1 else 0) +
          ∑ z, max 0 (P (z ∘ Equiv.swap K t) - P z)) :=
        sum_le_sum fun K _ => mul_le_mul_of_nonneg_left (hK K) (hq0 K)
    _ = ∑ z, P z * ∑ K, q K * (if 1 - α ≤ strangeWeight univ (S z) q K then 1 else 0) +
          ∑ K, q K * ∑ z, max 0 (P (z ∘ Equiv.swap K t) - P z) := by
        simp only [mul_add, sum_add_distrib]
        congr 1
        simp only [mul_sum]
        rw [sum_comm]
        exact sum_congr rfl fun z _ => sum_congr rfl fun K _ => by ring
    _ ≤ ∑ z, P z * (2 * α) + ∑ K, q K * ∑ z, max 0 (P (z ∘ Equiv.swap K t) - P z) :=
        add_le_add (sum_le_sum fun z _ => mul_le_mul_of_nonneg_left (hlem z) (hP0 z)) le_rfl
    _ = 2 * α * ∑ z, P z + ∑ K, q K * ∑ z, max 0 (P (z ∘ Equiv.swap K t) - P z) := by
        rw [mul_sum]
        congr 1
        exact sum_congr rfl fun z _ => by ring

/-- **Exchangeable data, weights fixed by position**: the swaps cost nothing, and the weighted
cross-fitted rule refuses with chance at most `2α`, whatever the fixed weights, so long as the test
position's is the largest. -/
theorem crossFit_weighted_exchangeable {Z β : Type*} [Fintype Z] [LinearOrder β]
    {S : (ι → Z) → ι → ι → β} (hS : PairEquivariant univ S) (P : (ι → Z) → ℝ)
    (hP0 : ∀ z, 0 ≤ P z) (q : ι → ℝ) (hq0 : ∀ i, 0 ≤ q i) (hq1 : ∑ i, q i = 1) (t : ι)
    (hqt : ∀ i, q i ≤ q t) (hPx : ∀ K z, P (z ∘ Equiv.swap K t) = P z) {α : ℝ} (hα : 0 ≤ α) :
    ∑ z, P z * (if 1 - α ≤ strangeWeight univ (S z) q t then 1 else 0) ≤ 2 * α * ∑ z, P z := by
  have h := crossFit_weighted_swap_tv hS P hP0 q hq0 hq1 t hqt hα
  simp only [hPx, sub_self, max_self, sum_const_zero, mul_zero, add_zero] at h
  exact h

/-- **The weighted jackknife+ refuses** the test unit when the upper ends `μ̂₋ᵢ(Xₜ) + Rᵢ` strictly
below `Yₜ`, or the lower ends strictly above it, carry at least `1 − α` of the total weight, the
test unit's own included. That is `Yₜ` outside Barber et al.'s weighted interval, with the test
unit's mass at `±∞`. -/
def JackRefusesW {Z X : Type*} (cov : Z → X) (resp : Z → ℝ) (α : ℝ) (C : Finset ι) (q : ι → ℝ)
    (alg : Finset ι → (ι → Z) → X → ℝ) (z : ι → Z) (t : ι) : Prop :=
  (1 - α) * ∑ j ∈ C, q j ≤ ∑ i ∈ (C.erase t).filter (fun i => alg (univ \ {i, t}) z (cov (z t)) +
      |resp (z i) - alg (univ \ {i, t}) z (cov (z i))| < resp (z t)), q i ∨
    (1 - α) * ∑ j ∈ C, q j ≤ ∑ i ∈ (C.erase t).filter (fun i => resp (z t) <
      alg (univ \ {i, t}) z (cov (z t)) - |resp (z i) - alg (univ \ {i, t}) z (cov (z i))|), q i

noncomputable instance {Z X : Type*} (cov : Z → X) (resp : Z → ℝ) (α : ℝ) (C : Finset ι)
    (q : ι → ℝ) (alg : Finset ι → (ι → Z) → X → ℝ) (z : ι → Z) (t : ι) :
    Decidable (JackRefusesW cov resp α C q alg z t) :=
  inferInstanceAs (Decidable (_ ∨ _))

/-- **The weighted jackknife+ refuses only a strange test unit**: the units its interval counts are
units the test unit is stranger than. -/
theorem strangeWeight_of_jackRefusesW {Z X : Type*} {cov : Z → X} {resp : Z → ℝ} {α : ℝ}
    {C : Finset ι} {q : ι → ℝ} (hq : ∀ j ∈ C, 0 ≤ q j) {alg : Finset ι → (ι → Z) → X → ℝ}
    {z : ι → Z} {t : ι} (h : JackRefusesW cov resp α C q alg z t) :
    (1 - α) * ∑ j ∈ C, q j ≤ strangeWeight C (pairScore (jackFit cov resp alg) z) q t := by
  have hpair : ∀ i, (univ \ {t, i} : Finset ι) = univ \ {i, t} := fun i => by rw [pair_comm]
  have hstr : ∀ j ∈ strangerThan C (pairScore (jackFit cov resp alg) z) t, 0 ≤ q j :=
    fun j hj => hq j (mem_filter.1 hj).1
  rcases h with h | h
  · refine h.trans (sum_le_sum_of_subset_of_nonneg (fun i hi => ?_) fun j hj _ => hstr j hj)
    rw [mem_filter, mem_erase] at hi
    rw [strangerThan, mem_filter]
    refine ⟨hi.1.2, ?_⟩
    simp only [pairScore, jackFit, hpair]
    rw [neg_lt_neg_iff]
    have h3 := le_abs_self (resp (z t) - alg (univ \ {i, t}) z (cov (z t)))
    linarith [hi.2]
  · refine h.trans (sum_le_sum_of_subset_of_nonneg (fun i hi => ?_) fun j hj _ => hstr j hj)
    rw [mem_filter, mem_erase] at hi
    rw [strangerThan, mem_filter]
    refine ⟨hi.1.2, ?_⟩
    simp only [pairScore, jackFit, hpair]
    rw [neg_lt_neg_iff]
    have h3 := le_abs_self (alg (univ \ {i, t}) z (cov (z t)) - resp (z t))
    rw [abs_sub_comm] at h3
    linarith [hi.2]

/-- **The nonexchangeable jackknife+** (Barber, Candès, Ramdas and Tibshirani 2023, Theorem 5), on
a finite law, with a symmetric regression fit: it misses with chance at most `2α` plus the
weighted swap distances `Σ_K q_K · TV(P, P ∘ swap(K, t))`. -/
theorem jackknifePlus_weighted {Z X : Type*} [Fintype Z] (cov : Z → X) (resp : Z → ℝ)
    {alg : Finset ι → (ι → Z) → X → ℝ} (halg : FitSymmetric univ alg) (P : (ι → Z) → ℝ)
    (hP0 : ∀ z, 0 ≤ P z) (q : ι → ℝ) (hq0 : ∀ i, 0 ≤ q i) (hq1 : ∑ i, q i = 1) (t : ι)
    (hqt : ∀ i, q i ≤ q t) {α : ℝ} (hα : 0 ≤ α) :
    ∑ z, P z * (if JackRefusesW cov resp α univ q alg z t then 1 else 0) ≤
      2 * α * ∑ z, P z + ∑ K, q K * ∑ z, max 0 (P (z ∘ Equiv.swap K t) - P z) := by
  classical
  have hsym : FitSymmetric univ (jackFit cov resp alg) := fun σ hσ T z => by
    funext p
    simp only [jackFit, halg σ hσ T z]
  refine le_trans (sum_le_sum fun z _ => mul_le_mul_of_nonneg_left ?_ (hP0 z))
    (crossFit_weighted_swap_tv (pairScore_equivariant hsym) P hP0 q hq0 hq1 t hqt hα)
  by_cases h : JackRefusesW cov resp α univ q alg z t
  · have := strangeWeight_of_jackRefusesW (fun j _ => hq0 j) h
    rw [hq1, mul_one] at this
    rw [if_pos h, if_pos this]
  · rw [if_neg h]
    split_ifs <;> norm_num

/-- **The tournament correction under weights**: the cross-fitted rule on the tournament's pair
score, with weights fixed by position, misses with chance at most `2α` plus the weighted swap
distances. Seeds that travel with the units are part of each unit's data. -/
theorem tournament_weighted {Z β α' : Type*} [Fintype Z] [LinearOrder α']
    {fit : Finset ι → (ι → Z) → β} (hfit : FitSymmetric univ fit) (g : β → Z → Multiset Z → α')
    (P : (ι → Z) → ℝ) (hP0 : ∀ z, 0 ≤ P z) (q : ι → ℝ) (hq0 : ∀ i, 0 ≤ q i)
    (hq1 : ∑ i, q i = 1) (t : ι) (hqt : ∀ i, q i ≤ q t) {α : ℝ} (hα : 0 ≤ α) :
    ∑ z, P z * (if 1 - α ≤ strangeWeight univ (pairScore (tournFit fit g) z) q t then 1 else 0) ≤
      2 * α * ∑ z, P z + ∑ K, q K * ∑ z, max 0 (P (z ∘ Equiv.swap K t) - P z) :=
  crossFit_weighted_swap_tv
    (pairScore_equivariant fun σ hσ T z => tournFit_comp g (hfit σ hσ T z)) P hP0 q hq0 hq1 t
    hqt hα

/-! #### Likelihood-ratio weights (JAWS) -/

/-- The weights travel with the units: reordering `C` reorders the comparisons and the weights. -/
theorem strangeWeight_comp {β : Type*} [LinearOrder β] {C : Finset ι} {σ : Equiv.Perm ι}
    (hσ : σ ∈ permsOn C) (s : ι → ι → β) (q : ι → ℝ) (u : ι) :
    strangeWeight C (fun i j => s (σ i) (σ j)) (q ∘ σ) u = strangeWeight C s q (σ u) := by
  unfold strangeWeight strangerThan
  exact sum_equiv σ (fun j => by simp only [mem_filter, mem_iff_of_mem_permsOn hσ])
    fun j _ => rfl

/-- **The cross-fitted rule with likelihood-ratio weights** (Prinster, Liu and Saria 2022, JAW's
Theorem 1, here for any pair-equivariant score). The weights `w(z_u)` travel with the units, the
test unit's own included; over the orderings of a data set the test unit lands on unit `u` with
weight `w(z_u)`, as under covariate shift (`w(x)`) or label shift (`w(y)`). Weighing each ordering
so, the weighted rule refuses at most `2α` of the weight: given the bag, the test unit is a strange
unit with chance their weight, at most `2α` by Lei's lemma. -/
theorem crossFit_likelihoodRatio {Z β : Type*} [LinearOrder β] {C : Finset ι} {t : ι}
    (ht : t ∈ C) {S : (ι → Z) → ι → ι → β} (hS : PairEquivariant C S) (z : ι → Z) (w : Z → ℝ)
    (hw : ∀ v, 0 ≤ w v) {α : ℝ} (hα : 0 ≤ α) :
    ∑ σ ∈ permsOn C, (if (1 - α) * ∑ j ∈ C, w (z (σ j)) ≤
        strangeWeight C (S (z ∘ σ)) (fun j => w (z (σ j))) t then w (z (σ t)) else 0) ≤
      2 * α * ∑ σ ∈ permsOn C, w (z (σ t)) := by
  have hc : (0 : ℝ) < C.card := by exact_mod_cast card_pos.2 ⟨t, ht⟩
  set f : ι → ℝ := fun u => if (1 - α) * ∑ j ∈ C, w (z j) ≤
    strangeWeight C (S z) (fun j => w (z j)) u then w (z u) else 0 with hf
  have hpt : ∀ σ ∈ permsOn C, (if (1 - α) * ∑ j ∈ C, w (z (σ j)) ≤
      strangeWeight C (S (z ∘ σ)) (fun j => w (z (σ j))) t then w (z (σ t)) else 0) =
      f (σ t) := by
    intro σ hσ
    have hW : ∑ j ∈ C, w (z (σ j)) = ∑ j ∈ C, w (z j) :=
      sum_equiv σ (fun j => (mem_iff_of_mem_permsOn hσ j).symm) fun j _ => rfl
    have hfun : S (z ∘ σ) = fun i j => S z (σ i) (σ j) :=
      funext fun i => funext fun j => hS σ hσ z i j
    rw [hW, hfun, show (fun j => w (z (σ j))) = (fun j => w (z j)) ∘ σ from rfl,
      strangeWeight_comp hσ]
  rw [sum_congr rfl hpt]
  have hlei := sum_weight_strange_le C (S z) (fun j => w (z j)) (fun j _ => hw (z j)) hα
  rw [sum_filter] at hlei
  refine (mul_le_mul_iff_of_pos_left hc).1 ?_
  rw [card_mul_sum_permsOn ht f, mul_left_comm (C.card : ℝ) (2 * α),
    card_mul_sum_permsOn ht (fun u => w (z u))]
  calc ((permsOn C).card : ℝ) * ∑ u ∈ C, f u ≤ (permsOn C).card * (2 * α * ∑ u ∈ C, w (z u)) :=
        mul_le_mul_of_nonneg_left hlei (Nat.cast_nonneg _)
    _ = 2 * α * ((permsOn C).card * ∑ u ∈ C, w (z u)) := by ring

/-- **JAW, the jackknife+ weighted by likelihood ratios** (Prinster, Liu and Saria 2022, Theorem 1).
With a symmetric fit, the interval with the units' weights `w(z_u)`, the test unit's own included
and its mass at `±∞`, misses with at most `2α` of the weight over the orderings. -/
theorem jackknifePlus_likelihoodRatio {Z X : Type*} {C : Finset ι} {t : ι} (ht : t ∈ C)
    (cov : Z → X) (resp : Z → ℝ) {alg : Finset ι → (ι → Z) → X → ℝ} (halg : FitSymmetric C alg)
    (z : ι → Z) (w : Z → ℝ) (hw : ∀ v, 0 ≤ w v) {α : ℝ} (hα : 0 ≤ α) :
    ∑ σ ∈ permsOn C, (if JackRefusesW cov resp α C (fun j => w (z (σ j))) alg (z ∘ σ) t
        then w (z (σ t)) else 0) ≤ 2 * α * ∑ σ ∈ permsOn C, w (z (σ t)) := by
  have hsym : FitSymmetric C (jackFit cov resp alg) := fun σ hσ T z => by
    funext p
    simp only [jackFit, halg σ hσ T z]
  refine le_trans (sum_le_sum fun σ _ => ?_)
    (crossFit_likelihoodRatio ht (pairScore_equivariant hsym) z w hw hα)
  by_cases h : JackRefusesW cov resp α C (fun j => w (z (σ j))) alg (z ∘ σ) t
  · rw [if_pos h, if_pos (strangeWeight_of_jackRefusesW (fun j _ => hw _) h)]
  · rw [if_neg h]
    split_ifs <;> linarith [hw (z (σ t))]

/-- **The tournament correction with likelihood-ratio weights**: at most `2α` of the weight over
the orderings, with a symmetric fit. -/
theorem tournament_likelihoodRatio {Z β α' : Type*} [LinearOrder α'] {C : Finset ι} {t : ι}
    (ht : t ∈ C) {fit : Finset ι → (ι → Z) → β} (hfit : FitSymmetric C fit)
    (g : β → Z → Multiset Z → α') (z : ι → Z) (w : Z → ℝ) (hw : ∀ v, 0 ≤ w v) {α : ℝ}
    (hα : 0 ≤ α) :
    ∑ σ ∈ permsOn C, (if (1 - α) * ∑ j ∈ C, w (z (σ j)) ≤
        strangeWeight C (pairScore (tournFit fit g) (z ∘ σ)) (fun j => w (z (σ j))) t
        then w (z (σ t)) else 0) ≤
      2 * α * ∑ σ ∈ permsOn C, w (z (σ t)) :=
  crossFit_likelihoodRatio ht
    (pairScore_equivariant fun σ hσ T z => tournFit_comp g (hfit σ hσ T z)) z w hw hα

/-! #### K-fold CV+ under weights -/

omit [Fintype ι] [DecidableEq ι] in
/-- **A weighted CV+ refusal is a refusal in each fold.** `U` holds the training units and the test
unit `t`; a unit is `beaten` when the test unit is stranger than it under its fold's fit. In a block
`C` of `U` scored by one fit `s` (a fold with the test unit), being unbeaten is scoring at most the
test unit. If the beaten units carry at least `1 − α` of `U`'s weight, the block's units scoring at
most the test unit, itself included, carry at most `α` of it. -/
theorem fold_refused_of_refused {U C : Finset ι} (hCU : C ⊆ U) (w : ι → ℝ)
    (hw : ∀ j ∈ U, 0 ≤ w j) (beaten : ι → Prop) [DecidablePred beaten] (s : ι → ℝ) (t : ι)
    (hb : ∀ j ∈ C, ¬ beaten j ↔ s j ≤ s t) {α : ℝ}
    (href : (1 - α) * ∑ j ∈ U, w j ≤ ∑ j ∈ U.filter beaten, w j) :
    ∑ j ∈ C.filter (fun j => s j ≤ s t), w j ≤ α * ∑ j ∈ U, w j := by
  have h1 : ∑ j ∈ C.filter (fun j => s j ≤ s t), w j ≤
      ∑ j ∈ U.filter (fun j => ¬ beaten j), w j := by
    refine sum_le_sum_of_subset_of_nonneg (fun j hj => ?_) fun j hj _ => hw j (mem_filter.1 hj).1
    rw [mem_filter] at hj ⊢
    exact ⟨hCU hj.1, (hb j hj.1).2 hj.2⟩
  have h2 := sum_filter_add_sum_filter_not U beaten w
  linarith

/-- **K-fold CV+ with likelihood-ratio weights, by one fold.** The block `C` is a fold with the test
unit, scored by that fold's fit `s`, which reads neither; `Wout ≥ 0` is the other folds' weight,
fixed given them. Over the reorderings of the block, each weighed by the test unit's `w`, the
block's refusal at level `α` of the whole weight, which every weighted CV+ refusal implies
(`fold_refused_of_refused`), carries at most `min(1, α (Wout + W_C) / W_C)` of the weight: the
weighted split rule in the block at that level. With the fold of largest weight this is the
bound to report; at equal weights it is slice 9's `ℓ/(m + 1)`. -/
theorem cvPlus_oneFold_likelihoodRatio {Z : Type*} {C : Finset ι} {t : ι} (ht : t ∈ C)
    (s : Z → ℝ) (z : ι → Z) (w : Z → ℝ) (hw : ∀ v, 0 ≤ w v) (Wout : ℝ) (hWout : 0 ≤ Wout)
    {α : ℝ} (hα : 0 ≤ α) :
    ∑ σ ∈ permsOn C, (if ∑ j ∈ C.filter (fun j => s (z (σ j)) ≤ s (z (σ t))), w (z (σ j)) ≤
        α * (Wout + ∑ j ∈ C, w (z j)) then w (z (σ t)) else 0) ≤
      min 1 (α * (Wout + ∑ j ∈ C, w (z j)) / ∑ j ∈ C, w (z j)) *
        ∑ σ ∈ permsOn C, w (z (σ t)) := by
  set WC := ∑ j ∈ C, w (z j) with hWC
  have hden : 0 ≤ ∑ σ ∈ permsOn C, w (z (σ t)) := sum_nonneg fun σ _ => hw _
  have hle1 : ∑ σ ∈ permsOn C, (if ∑ j ∈ C.filter (fun j => s (z (σ j)) ≤ s (z (σ t))),
      w (z (σ j)) ≤ α * (Wout + WC) then w (z (σ t)) else 0) ≤
      ∑ σ ∈ permsOn C, w (z (σ t)) :=
    sum_le_sum fun σ _ => by split_ifs <;> linarith [hw (z (σ t))]
  rcases (sum_nonneg fun i (_ : i ∈ C) => hw (z i)).eq_or_lt with h0 | hpos
  · have hz : ∀ i ∈ C, w (z i) = 0 :=
      (sum_eq_zero_iff_of_nonneg fun i _ => hw (z i)).1 h0.symm
    have hden0 : ∑ σ ∈ permsOn C, w (z (σ t)) = 0 :=
      sum_eq_zero fun σ hσ => hz _ ((mem_iff_of_mem_permsOn hσ t).2 ht)
    rw [hden0, mul_zero]
    rw [hden0] at hle1
    exact hle1
  · have hα' : 0 ≤ α * (Wout + WC) / WC :=
      div_nonneg (mul_nonneg hα (add_nonneg hWout hpos.le)) hpos.le
    have hcov := coverage_weighted (S := fun z => fun i => s (z i)) ht (fun σ _ z i => rfl) z w
      hw hα'
    have hpt : ∀ σ ∈ permsOn C, (if ∑ j ∈ C.filter (fun j => s (z (σ j)) ≤ s (z (σ t))),
        w (z (σ j)) ≤ α * (Wout + WC) then w (z (σ t)) else 0) =
        (if α * (Wout + WC) / WC < weightLe C (fun i => s ((z ∘ σ) i))
          (fun j => w (z (σ j)) / ∑ i ∈ C, w (z i)) t then 0 else w (z (σ t))) := by
      intro σ _
      have e : weightLe C (fun i => s ((z ∘ σ) i)) (fun j => w (z (σ j)) / ∑ i ∈ C, w (z i)) t =
          (∑ j ∈ C.filter (fun j => s (z (σ j)) ≤ s (z (σ t))), w (z (σ j))) / WC := by
        unfold weightLe
        rw [sum_div]
        rfl
      rw [e]
      by_cases h : ∑ j ∈ C.filter (fun j => s (z (σ j)) ≤ s (z (σ t))), w (z (σ j)) ≤
          α * (Wout + WC)
      · rw [if_pos h, if_neg (not_lt.2 ((div_le_div_iff_of_pos_right hpos).2 h))]
      · rw [if_neg h, if_pos ((div_lt_div_iff_of_pos_right hpos).2 (not_le.1 h))]
    rw [sum_congr rfl hpt]
    rw [sum_congr rfl hpt] at hle1
    rcases le_total 1 (α * (Wout + WC) / WC) with h1 | h1
    · rw [min_eq_left h1, one_mul]
      exact hle1
    · rw [min_eq_right h1]
      exact hcov

end WeightedCrossFit

/-! ### Slice 14's follow-up 3: CB-3's margin and inflated rule under weights

CB-3's statements with the weighted rule's `weightLe` in place of `rankLe`: the weighted fragile
units, the weighted inflated rule, and the band read from any ordering's estimates. -/

section WeightedMargin

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

omit [Fintype ι] in
/-- The weight of the units of `C` at least `δ j` below `s u`, `u` itself always counted. -/
noncomputable def weightShift (C : Finset ι) (s q δ : ι → ℝ) (u : ι) : ℝ :=
  ∑ j ∈ C.filter (fun j => j = u ∨ s j + δ j ≤ s u), q j

omit [Fintype ι] in
/-- **The weighted fragile units** under per-unit error bounds `ε`: accepted by the weighted rule on
the exact score, but with at most `α` of the weight at least `ε j + ε u` below, itself counted, so
estimates within their bounds can refuse it. -/
def WFragile (α : ℝ) (C : Finset ι) (s q ε : ι → ℝ) (u : ι) : Prop :=
  α < weightLe C s q u ∧ weightShift C s q (fun j => ε j + ε u) u ≤ α

omit [Fintype ι] [DecidableEq ι] in
/-- **The weighted inflated rule**: accept `u` when more than `α` of the weight has an estimate,
lowered by its bound, at most `u`'s raised by its own. -/
def WAcceptsInflated (α : ℝ) (C : Finset ι) (s q ε : ι → ℝ) (u : ι) : Prop :=
  α < ∑ j ∈ C.filter (fun j => s j - ε j ≤ s u + ε u), q j

omit [Fintype ι] in
/-- The weighted fragile band read from estimated scores and their bounds alone. -/
def WFragileEstimated (α : ℝ) (C : Finset ι) (s q ε : ι → ℝ) (u : ι) : Prop :=
  WAcceptsInflated α C s q ε u ∧ weightShift C s q (fun j => 2 * ε j + 2 * ε u) u ≤ α

omit [Fintype ι] in
noncomputable instance (α : ℝ) (C : Finset ι) (s q ε : ι → ℝ) (u : ι) :
    Decidable (WFragile α C s q ε u) :=
  inferInstanceAs (Decidable (_ ∧ _))

omit [Fintype ι] [DecidableEq ι] in
noncomputable instance (α : ℝ) (C : Finset ι) (s q ε : ι → ℝ) (u : ι) :
    Decidable (WAcceptsInflated α C s q ε u) :=
  inferInstanceAs (Decidable (_ < _))

omit [Fintype ι] in
noncomputable instance (α : ℝ) (C : Finset ι) (s q ε : ι → ℝ) (u : ι) :
    Decidable (WFragileEstimated α C s q ε u) :=
  inferInstanceAs (Decidable (_ ∧ _))

omit [Fintype ι] in
/-- Estimates within their bounds give each unit at least the exact score's shifted weight. -/
theorem weightShift_le_weightLe {C : Finset ι} {s ŝ ε q : ι → ℝ} (hq : ∀ j ∈ C, 0 ≤ q j)
    (hclose : ∀ j ∈ C, |ŝ j - s j| ≤ ε j) {u : ι} (hu : u ∈ C) :
    weightShift C s q (fun j => ε j + ε u) u ≤ weightLe C ŝ q u := by
  unfold weightShift weightLe
  refine sum_le_sum_of_subset_of_nonneg (fun j hj => ?_) fun j hj _ => hq j (mem_filter.1 hj).1
  rw [mem_filter] at hj ⊢
  refine ⟨hj.1, ?_⟩
  rcases hj.2 with rfl | h
  · exact le_rfl
  · have hj' := abs_le.1 (hclose j hj.1)
    have hu' := abs_le.1 (hclose u hu)
    linarith

omit [Fintype ι] in
/-- **CB-3 under weights, the deterministic part.** Estimates within their bounds refuse a unit the
exact score's weighted rule accepts only when that unit is weighted-fragile. -/
theorem waccepts_or_fragile {α : ℝ} {C : Finset ι} {s ŝ ε q : ι → ℝ} (hq : ∀ j ∈ C, 0 ≤ q j)
    (hclose : ∀ j ∈ C, |ŝ j - s j| ≤ ε j) {u : ι} (hu : u ∈ C) (h : α < weightLe C s q u) :
    α < weightLe C ŝ q u ∨ WFragile α C s q ε u := by
  by_cases hf : weightShift C s q (fun j => ε j + ε u) u ≤ α
  · exact Or.inr ⟨h, hf⟩
  · exact Or.inl (lt_of_lt_of_le (not_le.1 hf) (weightShift_le_weightLe hq hclose hu))

omit [Fintype ι] [DecidableEq ι] in
/-- The weighted inflated rule keeps every unit the exact score's weighted rule accepts. -/
theorem wAcceptsInflated_of_accepts {α : ℝ} {C : Finset ι} {s ŝ ε q : ι → ℝ}
    (hq : ∀ j ∈ C, 0 ≤ q j) (hclose : ∀ j ∈ C, |ŝ j - s j| ≤ ε j) {u : ι} (hu : u ∈ C)
    (h : α < weightLe C s q u) : WAcceptsInflated α C ŝ q ε u := by
  unfold WAcceptsInflated
  unfold weightLe at h
  refine lt_of_lt_of_le h (sum_le_sum_of_subset_of_nonneg (fun j hj => ?_)
    fun j hj _ => hq j (mem_filter.1 hj).1)
  rw [mem_filter] at hj ⊢
  have hj' := abs_le.1 (hclose j hj.1)
  have hu' := abs_le.1 (hclose u hu)
  exact ⟨hj.1, by linarith [hj.2]⟩

omit [Fintype ι] in
/-- **The engine's reading under weights.** A weighted-fragile unit of the exact score sits in the
band read from the estimated scores and their bounds alone. -/
theorem wFragileEstimated_of_fragile {α : ℝ} {C : Finset ι} {s ŝ ε q : ι → ℝ}
    (hq : ∀ j ∈ C, 0 ≤ q j) (hclose : ∀ j ∈ C, |ŝ j - s j| ≤ ε j) {u : ι} (hu : u ∈ C)
    (h : WFragile α C s q ε u) : WFragileEstimated α C ŝ q ε u := by
  refine ⟨wAcceptsInflated_of_accepts hq hclose hu h.1, le_trans ?_ h.2⟩
  unfold weightShift
  refine sum_le_sum_of_subset_of_nonneg (fun j hj => ?_) fun j hj _ => hq j (mem_filter.1 hj).1
  rw [mem_filter] at hj ⊢
  refine ⟨hj.1, ?_⟩
  rcases hj.2 with h' | h'
  · exact Or.inl h'
  · have hj' := abs_le.1 (hclose j hj.1)
    have hu' := abs_le.1 (hclose u hu)
    exact Or.inr (by linarith)

theorem weightShift_comp {C : Finset ι} {σ : Equiv.Perm ι} (hσ : σ ∈ permsOn C)
    (s q ε : ι → ℝ) (u : ι) :
    weightShift C (s ∘ σ) (q ∘ σ) (fun j => (ε ∘ σ) j + (ε ∘ σ) u) u =
      weightShift C s q (fun j => ε j + ε (σ u)) (σ u) := by
  unfold weightShift
  exact sum_equiv σ (fun j => by simp only [mem_filter, Function.comp_apply,
    mem_iff_of_mem_permsOn hσ, σ.injective.eq_iff]) fun j _ => rfl

theorem wFragile_comp {α : ℝ} {C : Finset ι} {σ : Equiv.Perm ι} (hσ : σ ∈ permsOn C)
    (s q ε : ι → ℝ) (u : ι) :
    WFragile α C (s ∘ σ) (q ∘ σ) (ε ∘ σ) u ↔ WFragile α C s q ε (σ u) := by
  simp only [WFragile, weightLe_comp hσ, weightShift_comp hσ]

/-- **CB-3 under weights, over orderings.** The exact score is equivariant, the weights `w(z_u)` and
the error bounds travel with the units, and the estimate `ŝ σ` for the ordering `σ` is anything
within those bounds. Weighing each ordering by the test unit's `w`, the estimate's weighted rule
refuses at most `α` plus `M`, the exact score's weighted-fragile units' share of the weight. -/
theorem coverage_close_weighted {Z : Type*} {C : Finset ι} {t : ι} (ht : t ∈ C)
    {S : (ι → Z) → ι → ℝ} (hS : Equivariant C S) (z : ι → Z) (w : Z → ℝ) (hw : ∀ v, 0 ≤ w v)
    (ŝ : Equiv.Perm ι → ι → ℝ) (ε : ι → ℝ)
    (hclose : ∀ σ ∈ permsOn C, ∀ j ∈ C, |ŝ σ j - S (z ∘ σ) j| ≤ ε (σ j)) {α : ℝ} (hα : 0 ≤ α) :
    ∑ σ ∈ permsOn C, (if α < weightLe C (ŝ σ) (fun j => w (z (σ j)) / ∑ i ∈ C, w (z i)) t then 0
        else w (z (σ t))) ≤
      (α + ∑ u ∈ C.filter (fun u => WFragile α C (S z) (fun j => w (z j) / ∑ i ∈ C, w (z i)) ε u),
        w (z u) / ∑ i ∈ C, w (z i)) * ∑ σ ∈ permsOn C, w (z (σ t)) := by
  set q : ι → ℝ := fun j => w (z j) / ∑ i ∈ C, w (z i) with hq
  have hq0 : ∀ j, 0 ≤ q j := fun j => div_nonneg (hw _) (sum_nonneg fun i _ => hw (z i))
  -- each ordering: a refusal by the estimate is a refusal by the exact score or a fragile test unit
  have hpt : ∀ σ ∈ permsOn C, (if α < weightLe C (ŝ σ) (fun j => w (z (σ j)) / ∑ i ∈ C, w (z i)) t
      then 0 else w (z (σ t))) ≤
      (if α < weightLe C (S (z ∘ σ)) (fun j => w (z (σ j)) / ∑ i ∈ C, w (z i)) t then 0
        else w (z (σ t))) + (fun u => if WFragile α C (S z) q ε u then w (z u) else 0) (σ t) := by
    intro σ hσ
    have hfun : S (z ∘ σ) = S z ∘ σ := funext fun i => hS σ hσ z i
    have hqσ : (fun j => w (z (σ j)) / ∑ i ∈ C, w (z i)) = q ∘ σ := rfl
    simp only
    by_cases he : α < weightLe C (ŝ σ) (fun j => w (z (σ j)) / ∑ i ∈ C, w (z i)) t
    · rw [if_pos he]
      split_ifs <;> linarith [hw (z (σ t))]
    · rw [if_neg he]
      by_cases hx : α < weightLe C (S (z ∘ σ)) (fun j => w (z (σ j)) / ∑ i ∈ C, w (z i)) t
      · rw [if_pos hx, zero_add]
        rcases waccepts_or_fragile (ε := ε ∘ σ) (fun j _ => hq0 (σ j)) (hclose σ hσ) ht hx with
          h | h
        · exact absurd h he
        · rw [hfun, hqσ, wFragile_comp hσ] at h
          rw [if_pos h]
      · rw [if_neg hx]
        split_ifs <;> linarith [hw (z (σ t))]
  have hsplit := sum_le_sum hpt
  rw [sum_add_distrib] at hsplit
  have hexact := coverage_weighted ht hS z w hw hα
  have hc : (0 : ℝ) < C.card := by exact_mod_cast card_pos.2 ⟨t, ht⟩
  -- the fragile term over orderings is the fragile units' share of the weight
  have hfrag : ∑ σ ∈ permsOn C, (fun u => if WFragile α C (S z) q ε u then w (z u) else 0) (σ t) =
      (∑ u ∈ C.filter (fun u => WFragile α C (S z) q ε u), q u) *
        ∑ σ ∈ permsOn C, w (z (σ t)) := by
    refine mul_left_cancel₀ hc.ne' ?_
    rw [card_mul_sum_permsOn ht (fun u => if WFragile α C (S z) q ε u then w (z u) else 0),
      mul_left_comm (C.card : ℝ), card_mul_sum_permsOn ht (fun u => w (z u)), ← sum_filter]
    rcases (sum_nonneg fun i (_ : i ∈ C) => hw (z i)).eq_or_lt with h0 | hpos
    · have hz : ∀ i ∈ C, w (z i) = 0 :=
        (sum_eq_zero_iff_of_nonneg fun i _ => hw (z i)).1 h0.symm
      rw [sum_eq_zero fun u hu => hz u (mem_filter.1 hu).1, ← h0]
      ring
    · rw [hq, ← sum_div]
      field_simp
  rw [hfrag] at hsplit
  linarith

/-- **The weighted inflated rule over orderings** keeps the weighted rule's guarantee, at the price
of a wider set. -/
theorem coverage_inflated_weighted {Z : Type*} {C : Finset ι} {t : ι} (ht : t ∈ C)
    {S : (ι → Z) → ι → ℝ} (hS : Equivariant C S) (z : ι → Z) (w : Z → ℝ) (hw : ∀ v, 0 ≤ w v)
    (ŝ ε : Equiv.Perm ι → ι → ℝ)
    (hclose : ∀ σ ∈ permsOn C, ∀ j ∈ C, |ŝ σ j - S (z ∘ σ) j| ≤ ε σ j) {α : ℝ} (hα : 0 ≤ α) :
    ∑ σ ∈ permsOn C, (if WAcceptsInflated α C (ŝ σ)
        (fun j => w (z (σ j)) / ∑ i ∈ C, w (z i)) (ε σ) t then 0 else w (z (σ t))) ≤
      α * ∑ σ ∈ permsOn C, w (z (σ t)) := by
  have hq0 : ∀ j, 0 ≤ w (z j) / ∑ i ∈ C, w (z i) :=
    fun j => div_nonneg (hw _) (sum_nonneg fun i _ => hw (z i))
  refine le_trans (sum_le_sum fun σ hσ => ?_) (coverage_weighted ht hS z w hw hα)
  by_cases hx : α < weightLe C (S (z ∘ σ)) (fun j => w (z (σ j)) / ∑ i ∈ C, w (z i)) t
  · rw [if_pos hx, if_pos (wAcceptsInflated_of_accepts (fun j _ => hq0 (σ j)) (hclose σ hσ) ht hx)]
  · rw [if_neg hx]
    split_ifs <;> linarith [hw (z (σ t))]

/-- **The weighted margin read from the estimate.** For every ordering, the exact score's
weighted-fragile units weigh at most the band read from that ordering's estimates, the weights and
bounds travelling with the units, so the engine can report `α + M̂` with `M̂` read from what it
holds. -/
theorem sum_wfragile_le_estimated {Z : Type*} {C : Finset ι} {S : (ι → Z) → ι → ℝ}
    (hS : Equivariant C S) (z : ι → Z) (q : ι → ℝ) (hq : ∀ j ∈ C, 0 ≤ q j) {σ : Equiv.Perm ι}
    (hσ : σ ∈ permsOn C) (ŝ ε : ι → ℝ) (hclose : ∀ j ∈ C, |ŝ j - S (z ∘ σ) j| ≤ ε (σ j))
    (α : ℝ) :
    ∑ u ∈ C.filter (fun u => WFragile α C (S z) q ε u), q u ≤
      ∑ u ∈ C.filter (fun u => WFragileEstimated α C ŝ (q ∘ σ) (ε ∘ σ) u), q (σ u) := by
  have hfun : S (z ∘ σ) = S z ∘ σ := funext fun i => hS σ hσ z i
  have hqσ : ∀ j ∈ C, 0 ≤ (q ∘ σ) j := fun j hj => hq (σ j) ((mem_iff_of_mem_permsOn hσ j).2 hj)
  have heq : ∑ u ∈ C.filter (fun u => WFragile α C (S z) q ε u), q u =
      ∑ u ∈ C.filter (fun u => WFragile α C (S (z ∘ σ)) (q ∘ σ) (ε ∘ σ) u), q (σ u) := by
    symm
    refine sum_equiv σ (fun u => ?_) fun u _ => rfl
    simp only [mem_filter, hfun, wFragile_comp hσ, mem_iff_of_mem_permsOn hσ]
  rw [heq]
  refine sum_le_sum_of_subset_of_nonneg (fun u hu => ?_) fun u hu _ => hqσ u (mem_filter.1 hu).1
  rw [mem_filter] at hu ⊢
  exact ⟨hu.1, wFragileEstimated_of_fragile hqσ hclose hu.1 hu.2⟩

end WeightedMargin

/-! ### Slice 14's follow-ups: K-fold CV+ under position weights, by merging the folds

The weighted counterpart of slice 9's merged bound. A weighted CV+ refusal makes a fixed-weight
average of the folds' weighted p-values small; averaging super-uniform p-values with fixed weights
costs a factor `2`, the arithmetic mean's factor in slice 9's merged bound, here with real-valued
p-values through a finite layer cake. -/

section WeightedMerge

/-- **A finite layer cake.** Masses `m ≥ 0` on finitely many levels `T`, each at most `b`, whose
cumulative mass up to every level is at most that level: `Σ m_v (b − v) ≤ b²/2`. -/
theorem sum_mass_mul_sub_le (T : Finset ℝ) (m : ℝ → ℝ) (hm : ∀ v ∈ T, 0 ≤ m v) (b : ℝ)
    (hb : ∀ v ∈ T, v ≤ b) (hcum : ∀ v ∈ T, ∑ w ∈ T.filter (· ≤ v), m w ≤ v) :
    ∑ v ∈ T, m v * (b - v) ≤ b ^ 2 / 2 := by
  induction T using Finset.induction_on_max generalizing b with
  | empty =>
    simp only [sum_empty]
    positivity
  | insert M T hM ih =>
    have hMT : M ∉ T := fun h => lt_irrefl M (hM M h)
    have hT : ∀ v ∈ T, (insert M T).filter (· ≤ v) = T.filter (· ≤ v) := fun v hv => by
      rw [filter_insert, if_neg (not_le.2 (hM v hv))]
    have ih' := ih (fun v hv => hm v (mem_insert_of_mem hv)) M (fun v hv => (hM v hv).le)
      (fun v hv => by rw [← hT v hv]; exact hcum v (mem_insert_of_mem hv))
    have hGM : ∑ w ∈ (insert M T).filter (· ≤ M), m w ≤ M := hcum M (mem_insert_self M T)
    have hfM : (insert M T).filter (· ≤ M) = insert M T := filter_true_of_mem fun v hv => by
      rcases mem_insert.1 hv with rfl | hv'
      · exact le_rfl
      · exact (hM v hv').le
    rw [hfM, sum_insert hMT] at hGM
    have hbM : M ≤ b := hb M (mem_insert_self M T)
    have hmM := hm M (mem_insert_self M T)
    have hsplit : ∑ v ∈ T, m v * (b - v) =
        ∑ v ∈ T, m v * (M - v) + (b - M) * ∑ v ∈ T, m v := by
      rw [mul_sum, ← sum_add_distrib]
      exact sum_congr rfl fun v _ => by ring
    rw [sum_insert hMT, hsplit]
    nlinarith [mul_le_mul_of_nonneg_left hGM (sub_nonneg.2 hbM), sq_nonneg (b - M)]

/-- **The layer cake for a super-uniform p-value**: `E(b − π)⁺ ≤ b²/2` when `π ≥ 0` and
`P(π ≤ x) ≤ x` at every level `0 ≤ x < b`. -/
theorem sum_mul_max_sub_le {Ω : Type*} [Fintype Ω] (P : Ω → ℝ) (hP : ∀ ω, 0 ≤ P ω) (π : Ω → ℝ)
    (hπ : ∀ ω, 0 ≤ π ω) (b : ℝ)
    (hsu : ∀ x, 0 ≤ x → x < b → ∑ ω, P ω * (if π ω ≤ x then 1 else 0) ≤ x) :
    ∑ ω, P ω * max 0 (b - π ω) ≤ b ^ 2 / 2 := by
  classical
  set s := univ.filter fun ω => π ω < b with hs
  set m : ℝ → ℝ := fun v => ∑ ω ∈ s.filter (fun ω => π ω = v), P ω with hmdef
  have h1 : ∑ ω, P ω * max 0 (b - π ω) = ∑ ω ∈ s, P ω * (b - π ω) := by
    rw [hs, sum_filter]
    refine sum_congr rfl fun ω _ => ?_
    split_ifs with h
    · rw [max_eq_right (by linarith)]
    · rw [max_eq_left (by linarith [not_lt.1 h]), mul_zero]
  have h2 : ∑ ω ∈ s, P ω * (b - π ω) = ∑ v ∈ s.image π, m v * (b - v) := by
    refine (sum_image' (fun ω => P ω * (b - π ω)) fun i _ => ?_).symm
    change (∑ ω ∈ s.filter (fun ω => π ω = π i), P ω) * (b - π i) =
      ∑ j ∈ s with π j = π i, P j * (b - π j)
    rw [sum_mul]
    exact sum_congr rfl fun j hj => by rw [(mem_filter.1 hj).2]
  rw [h1, h2]
  refine sum_mass_mul_sub_le _ m (fun v _ => sum_nonneg fun ω _ => hP ω) b
    (fun v hv => ?_) (fun v hv => ?_)
  · obtain ⟨ω, hω, rfl⟩ := mem_image.1 hv
    exact (mem_filter.1 hω).2.le
  · obtain ⟨ω₀, hω₀, rfl⟩ := mem_image.1 hv
    have hvb : π ω₀ < b := (mem_filter.1 hω₀).2
    have himg : (s.image π).filter (· ≤ π ω₀) = (s.filter fun ω => π ω ≤ π ω₀).image π := by
      ext v
      simp only [mem_filter, mem_image]
      constructor
      · rintro ⟨⟨ω, hω, rfl⟩, hle⟩
        exact ⟨ω, ⟨hω, hle⟩, rfl⟩
      · rintro ⟨ω, ⟨hω, hle⟩, rfl⟩
        exact ⟨⟨ω, hω, rfl⟩, hle⟩
    have hcum : ∑ w ∈ (s.image π).filter (· ≤ π ω₀), m w =
        ∑ ω ∈ s.filter (fun ω => π ω ≤ π ω₀), P ω := by
      rw [himg]
      refine sum_image' P fun i hi => ?_
      change ∑ ω ∈ s.filter (fun ω => π ω = π i), P ω =
        ∑ j ∈ s.filter (fun ω => π ω ≤ π ω₀) with π j = π i, P j
      have hi' : π i ≤ π ω₀ := (mem_filter.1 hi).2
      refine sum_congr ?_ fun _ _ => rfl
      ext j
      simp only [mem_filter]
      constructor
      · rintro ⟨hj, h⟩
        exact ⟨⟨hj, by rw [h]; exact hi'⟩, h⟩
      · rintro ⟨⟨hj, _⟩, h⟩
        exact ⟨hj, h⟩
    rw [hcum]
    have hall : ∑ ω ∈ s.filter (fun ω => π ω ≤ π ω₀), P ω =
        ∑ ω, P ω * (if π ω ≤ π ω₀ then 1 else 0) := by
      simp only [mul_ite, mul_one, mul_zero]
      rw [← sum_filter]
      refine sum_congr ?_ fun _ _ => rfl
      ext ω
      simp only [hs, mem_filter, mem_univ, true_and]
      constructor
      · rintro ⟨_, h⟩
        exact h
      · intro h
        exact ⟨lt_of_le_of_lt h hvb, h⟩
    rw [hall]
    exact hsu _ (hπ ω₀) hvb

/-- **Averaging super-uniform p-values with fixed weights costs a factor `2`** (the arithmetic
mean's factor, as in slice 9's merged bound): with `λ ≥ 0` summing to one and each `π_k ≥ 0`
super-uniform below `2c`, the average is at most `c` with chance at most `2c`. Where the average
is at most `c`, `Σ λ_k (2c − π_k)⁺ ≥ c`; each term averages at most `(2c)²/2` by the layer
cake. -/
theorem merge_le_two {Ω κ : Type*} [Fintype Ω] (P : Ω → ℝ) (hP : ∀ ω, 0 ≤ P ω) (K : Finset κ)
    (lam : κ → ℝ) (hlam : ∀ k ∈ K, 0 ≤ lam k) (hlam1 : ∑ k ∈ K, lam k = 1) (π : κ → Ω → ℝ)
    (hπ : ∀ k ∈ K, ∀ ω, 0 ≤ π k ω) {c : ℝ} (hc : 0 < c)
    (hsu : ∀ k ∈ K, ∀ x, 0 ≤ x → x < 2 * c → ∑ ω, P ω * (if π k ω ≤ x then 1 else 0) ≤ x) :
    ∑ ω, P ω * (if ∑ k ∈ K, lam k * π k ω ≤ c then 1 else 0) ≤ 2 * c := by
  have hpt : ∀ ω, c * (if ∑ k ∈ K, lam k * π k ω ≤ c then 1 else 0) ≤
      ∑ k ∈ K, lam k * max 0 (2 * c - π k ω) := by
    intro ω
    split_ifs with h
    · have : ∑ k ∈ K, lam k * (2 * c - π k ω) ≤ ∑ k ∈ K, lam k * max 0 (2 * c - π k ω) :=
        sum_le_sum fun k hk => mul_le_mul_of_nonneg_left (le_max_right _ _) (hlam k hk)
      have e : ∑ k ∈ K, lam k * (2 * c - π k ω) = 2 * c - ∑ k ∈ K, lam k * π k ω := by
        simp only [mul_sub, sum_sub_distrib]
        rw [← sum_mul, hlam1]
        ring
      linarith
    · rw [mul_zero]
      exact sum_nonneg fun k hk => mul_nonneg (hlam k hk) (le_max_left _ _)
  have hlc : ∀ k ∈ K, ∑ ω, P ω * max 0 (2 * c - π k ω) ≤ (2 * c) ^ 2 / 2 := fun k hk =>
    sum_mul_max_sub_le P hP (π k) (hπ k hk) (2 * c) (hsu k hk)
  have hsum : c * ∑ ω, P ω * (if ∑ k ∈ K, lam k * π k ω ≤ c then 1 else 0) ≤
      ∑ k ∈ K, lam k * ((2 * c) ^ 2 / 2) := by
    calc c * ∑ ω, P ω * (if ∑ k ∈ K, lam k * π k ω ≤ c then 1 else 0)
        = ∑ ω, P ω * (c * (if ∑ k ∈ K, lam k * π k ω ≤ c then 1 else 0)) := by
          rw [mul_sum]
          exact sum_congr rfl fun ω _ => by ring
      _ ≤ ∑ ω, P ω * ∑ k ∈ K, lam k * max 0 (2 * c - π k ω) :=
          sum_le_sum fun ω _ => mul_le_mul_of_nonneg_left (hpt ω) (hP ω)
      _ = ∑ k ∈ K, lam k * ∑ ω, P ω * max 0 (2 * c - π k ω) := by
          simp only [mul_sum]
          rw [sum_comm]
          exact sum_congr rfl fun k _ => sum_congr rfl fun ω _ => by ring
      _ ≤ ∑ k ∈ K, lam k * ((2 * c) ^ 2 / 2) :=
          sum_le_sum fun k hk => mul_le_mul_of_nonneg_left (hlc k hk) (hlam k hk)
  rw [← sum_mul, hlam1, one_mul] at hsum
  nlinarith

variable {ι κ : Type*}

/-- **A weighted CV+ refusal makes the folds' average small.** Training units `U` in folds
`fold : ι → κ` over `K`, weights `w ≥ 0` fixed by position and the test unit's `wt > 0`; a unit is
`beaten` when the test unit is stranger than it under its fold's fit. Fold `k`'s weighted p-value
is `π_k = (W_k⁻ + wt)/(W_k + wt)`, `W_k⁻` its unbeaten weight. If the beaten units carry at least
`1 − α` of `W + wt`, then `Σ_k λ_k π_k ≤ c` with `λ_k = (W_k + wt)/(W + |K| wt)` and
`c = (α (W + wt) + (|K| − 1) wt)/(W + |K| wt)`. -/
theorem cv_merged_of_refused (U : Finset ι) (K : Finset κ) (hK : K.Nonempty) (fold : ι → κ)
    (hfold : ∀ i ∈ U, fold i ∈ K) (w : ι → ℝ) (hw : ∀ i ∈ U, 0 ≤ w i) (wt : ℝ) (hwt : 0 < wt)
    (beaten : ι → Prop) [DecidablePred beaten] [DecidableEq κ] {α : ℝ}
    (href : (1 - α) * (∑ i ∈ U, w i + wt) ≤ ∑ i ∈ U.filter beaten, w i) :
    ∑ k ∈ K, (∑ i ∈ U.filter (fun i => fold i = k), w i + wt) / (∑ i ∈ U, w i + K.card * wt) *
        ((∑ i ∈ (U.filter fun i => fold i = k).filter (fun i => ¬ beaten i), w i + wt) /
          (∑ i ∈ U.filter (fun i => fold i = k), w i + wt)) ≤
      (α * (∑ i ∈ U, w i + wt) + (K.card - 1) * wt) / (∑ i ∈ U, w i + K.card * wt) := by
  have hW : 0 ≤ ∑ i ∈ U, w i := sum_nonneg hw
  have hD : 0 < ∑ i ∈ U, w i + K.card * wt := by
    have h1 : (1 : ℝ) ≤ K.card := by exact_mod_cast card_pos.2 hK
    nlinarith
  have hfk : ∀ k, 0 < ∑ i ∈ U.filter (fun i => fold i = k), w i + wt := fun k =>
    add_pos_of_nonneg_of_pos (sum_nonneg fun i hi => hw i (mem_filter.1 hi).1) hwt
  have hterm : ∀ k ∈ K, (∑ i ∈ U.filter (fun i => fold i = k), w i + wt) /
      (∑ i ∈ U, w i + K.card * wt) *
        ((∑ i ∈ (U.filter fun i => fold i = k).filter (fun i => ¬ beaten i), w i + wt) /
          (∑ i ∈ U.filter (fun i => fold i = k), w i + wt)) =
      (∑ i ∈ (U.filter fun i => fold i = k).filter (fun i => ¬ beaten i), w i + wt) /
        (∑ i ∈ U, w i + K.card * wt) := fun k _ => by
    have := (hfk k).ne'
    field_simp
  rw [sum_congr rfl hterm, ← sum_div]
  refine div_le_div_of_nonneg_right ?_ hD.le
  rw [sum_add_distrib, sum_const, nsmul_eq_mul]
  have hfib : ∑ k ∈ K, ∑ i ∈ (U.filter fun i => fold i = k).filter (fun i => ¬ beaten i), w i =
      ∑ i ∈ U.filter (fun i => ¬ beaten i), w i := by
    rw [← sum_fiberwise_of_maps_to (s := U.filter fun i => ¬ beaten i) (t := K) (g := fold)
      (fun i hi => hfold i (mem_filter.1 hi).1)]
    exact sum_congr rfl fun k _ =>
      sum_congr (by ext i; simp only [mem_filter]; tauto) fun _ _ => rfl
  rw [hfib]
  have hsplit := sum_filter_add_sum_filter_not U beaten w
  linarith

/-- **K-fold CV+ under position weights, by merging the folds** (our statement, the weighted form
of slice 9's merged bound). On a finite law `P` of the outcomes `ω`, the beaten units and so the
folds' weighted p-values depend on `ω`; the weights are fixed by position, the test unit's `wt > 0`.
With each fold's weighted p-value super-uniform below `2c` (on exchangeable data, Barber et al.'s
Theorem 2 in that fold, `refused_le_swap_tv`, with the test position's weight the largest in it),
CV+ refuses with chance at most `2c`, `c = (α (W + wt) + (|K| − 1) wt)/(W + |K| wt)`. At equal
weights `2c` is slice 9's `2(ℓ + K − 1)/(n + K)`. -/
theorem cvPlus_merged_weighted {Ω : Type*} [Fintype Ω] (P : Ω → ℝ) (hP : ∀ ω, 0 ≤ P ω)
    (U : Finset ι) (K : Finset κ) (hK : K.Nonempty) (fold : ι → κ) (hfold : ∀ i ∈ U, fold i ∈ K)
    (w : ι → ℝ) (hw : ∀ i ∈ U, 0 ≤ w i) (wt : ℝ) (hwt : 0 < wt) (beaten : Ω → ι → Prop)
    [∀ ω, DecidablePred (beaten ω)] [DecidableEq κ] {α : ℝ} (hα : 0 < α)
    (hsu : ∀ k ∈ K, ∀ x, 0 ≤ x → x < 2 * ((α * (∑ i ∈ U, w i + wt) + (K.card - 1) * wt) /
        (∑ i ∈ U, w i + K.card * wt)) →
      ∑ ω, P ω * (if (∑ i ∈ (U.filter fun i => fold i = k).filter (fun i => ¬ beaten ω i), w i +
          wt) / (∑ i ∈ U.filter (fun i => fold i = k), w i + wt) ≤ x then 1 else 0) ≤ x) :
    ∑ ω, P ω * (if (1 - α) * (∑ i ∈ U, w i + wt) ≤ ∑ i ∈ U.filter (beaten ω), w i then 1 else 0) ≤
      2 * ((α * (∑ i ∈ U, w i + wt) + (K.card - 1) * wt) / (∑ i ∈ U, w i + K.card * wt)) := by
  set D := ∑ i ∈ U, w i + K.card * wt with hDdef
  set c := (α * (∑ i ∈ U, w i + wt) + (K.card - 1) * wt) / D with hcdef
  have hW : 0 ≤ ∑ i ∈ U, w i := sum_nonneg hw
  have h1 : (1 : ℝ) ≤ K.card := by exact_mod_cast card_pos.2 hK
  have hD : 0 < D := by nlinarith
  have hc : 0 < c := by
    rw [hcdef]
    refine div_pos ?_ hD
    nlinarith
  set lam : κ → ℝ := fun k => (∑ i ∈ U.filter (fun i => fold i = k), w i + wt) / D with hlamdef
  have hlam : ∀ k ∈ K, 0 ≤ lam k := fun k _ =>
    div_nonneg (add_nonneg (sum_nonneg fun i hi => hw i (mem_filter.1 hi).1) hwt.le) hD.le
  have hlam1 : ∑ k ∈ K, lam k = 1 := by
    rw [hlamdef]
    simp only
    rw [← sum_div, sum_add_distrib, sum_const, nsmul_eq_mul,
      sum_fiberwise_of_maps_to (s := U) (t := K) (g := fold) hfold, div_self hD.ne']
  have hπ : ∀ k ∈ K, ∀ ω, 0 ≤ (∑ i ∈ (U.filter fun i => fold i = k).filter (fun i => ¬ beaten ω i),
      w i + wt) / (∑ i ∈ U.filter (fun i => fold i = k), w i + wt) := fun k _ ω =>
    div_nonneg
      (add_nonneg (sum_nonneg fun i hi => hw i (mem_filter.1 (mem_filter.1 hi).1).1) hwt.le)
      (add_nonneg (sum_nonneg fun i hi => hw i (mem_filter.1 hi).1) hwt.le)
  refine le_trans (sum_le_sum fun ω _ => mul_le_mul_of_nonneg_left ?_ (hP ω))
    (merge_le_two P hP K lam hlam hlam1 (fun k ω =>
      (∑ i ∈ (U.filter fun i => fold i = k).filter (fun i => ¬ beaten ω i), w i + wt) /
        (∑ i ∈ U.filter (fun i => fold i = k), w i + wt)) hπ hc hsu)
  by_cases h : (1 - α) * (∑ i ∈ U, w i + wt) ≤ ∑ i ∈ U.filter (beaten ω), w i
  · rw [if_pos h, if_pos (cv_merged_of_refused U K hK fold hfold w hw wt hwt (beaten ω) h)]
  · rw [if_neg h]
    split_ifs <;> norm_num

end WeightedMerge

/-! ### Slice 14's follow-ups: the tournament across folds under likelihood-ratio weights

Slice 9's CV-b with weights. Units play only across folds; a unit is refused when the units it beats
carry at least `1 − α` of the weight outside its fold plus its own. The refused units weigh at most
`2αW` plus the heaviest fold, and over the fold and slot permutations, each weighed by the test
unit's likelihood ratio, the test unit is refused with at most `2α + W_max/W` of the weight. -/

section WeightedAcross

variable {ι β : Type*} [DecidableEq β]

/-- The weight of the units of `C` in other folds than `u`'s that `u` is stranger than. -/
def strangeAcrossWeight {α' : Type*} [LinearOrder α'] (C : Finset ι) (f : ι → β)
    (s : ι → ι → α') (q : ι → ℝ) (u : ι) : ℝ :=
  ∑ j ∈ strangerAcross C f s u, q j

/-- **The weighted counting core across folds** (our statement; Barber et al.'s Step 1 for CV+,
weighted, as Lei's lemma weighs the jackknife+'s). Units play only across folds. `u` is refused when
the units it beats carry at least `1 − α` of the weight outside its fold plus its own. Then the
refused units weigh at most `2αW + W_max`, `W_max` bounding every fold's weight: each refused unit
loses at most `α W` across folds, and among the refused units each comparison across folds goes at
most one way. -/
theorem sum_weight_refused_across_le {α' : Type*} [LinearOrder α'] (C : Finset ι) (f : ι → β)
    (s : ι → ι → α') (q : ι → ℝ) (hq : ∀ j ∈ C, 0 ≤ q j) {α Wmax : ℝ} (hα0 : 0 ≤ α)
    (hα1 : α ≤ 1) (hW0 : 0 ≤ Wmax)
    (hWmax : ∀ u ∈ C, ∑ j ∈ C.filter (fun j => f j = f u), q j ≤ Wmax) :
    ∑ u ∈ C.filter (fun u => (1 - α) * (∑ j ∈ C, q j - ∑ j ∈ C.filter (fun j => f j = f u), q j +
        q u) ≤ strangeAcrossWeight C f s q u), q u ≤ 2 * α * ∑ j ∈ C, q j + Wmax := by
  classical
  set F := C.filter (fun u => (1 - α) * (∑ j ∈ C, q j -
    ∑ j ∈ C.filter (fun j => f j = f u), q j + q u) ≤ strangeAcrossWeight C f s q u) with hF
  set W := ∑ j ∈ C, q j with hW
  set Q := ∑ u ∈ F, q u with hQ
  have hFC : F ⊆ C := filter_subset _ _
  have hqF : ∀ j ∈ F, 0 ≤ q j := fun j hj => hq j (hFC hj)
  set a : ι → ℝ := fun u => ∑ j ∈ F, if f j ≠ f u ∧ s u j < s j u then q j else 0 with ha
  set Qf : ι → ℝ := fun u => ∑ j ∈ F, if f j = f u then q j else 0 with hQf
  set Wf : ι → ℝ := fun u => ∑ j ∈ C.filter (fun j => f j = f u), q j with hWf
  -- each refused unit: the units it beats in F, against the weight it could beat outside F
  have hunit : ∀ u ∈ F, Q - Qf u - α * W ≤ a u := by
    intro u hu
    have h1 := (mem_filter.1 hu).2
    have hrest : ∑ j ∈ C \ F, (if f j ≠ f u ∧ s u j < s j u then q j else 0) ≤
        ∑ j ∈ C \ F, (if f j ≠ f u then q j else 0) :=
      sum_le_sum fun j hj => by
        have := hq j (mem_sdiff.1 hj).1
        split_ifs with h h' <;> first | exact le_rfl | exact absurd h.1 h' | exact this
    have hsw : strangeAcrossWeight C f s q u =
        ∑ j ∈ C \ F, (if f j ≠ f u ∧ s u j < s j u then q j else 0) + a u := by
      unfold strangeAcrossWeight strangerAcross
      rw [sum_filter, ← sum_sdiff hFC]
    have hout : ∑ j ∈ C \ F, (if f j ≠ f u then q j else 0) = (W - Q) - (Wf u - Qf u) := by
      have e1 : ∑ j ∈ C \ F, (if f j ≠ f u then q j else 0) =
          ∑ j ∈ C \ F, q j - ∑ j ∈ C \ F, (if f j = f u then q j else 0) := by
        rw [← sum_sub_distrib]
        exact sum_congr rfl fun j _ => by
          by_cases h : f j = f u <;> simp [h]
      have e2 : ∑ j ∈ C \ F, q j = W - Q := by rw [sum_sdiff_eq_sub hFC]
      have e3 : ∑ j ∈ C \ F, (if f j = f u then q j else 0) = Wf u - Qf u := by
        have := sum_sdiff hFC (f := fun j => if f j = f u then q j else 0)
        have hWf' : Wf u = ∑ j ∈ C, (if f j = f u then q j else 0) := by
          rw [hWf]
          simp only
          rw [sum_filter]
        rw [hWf']
        change _ = _ - ∑ j ∈ F, (if f j = f u then q j else 0)
        linarith
      rw [e1, e2, e3]
    have hWfu : 0 ≤ Wf u := sum_nonneg fun j hj => hq j (mem_filter.1 hj).1
    have hqu : 0 ≤ q u := hqF u hu
    have hst := h1
    rw [hsw] at hst
    have h3 : (1 - α) * q u ≥ 0 := mul_nonneg (by linarith) hqu
    nlinarith [mul_nonneg hα0 hWfu]
  -- summed with the refused units' weights
  have hsum : ∑ u ∈ F, q u * (Q - Qf u - α * W) ≤ ∑ u ∈ F, q u * a u :=
    sum_le_sum fun u hu => mul_le_mul_of_nonneg_left (hunit u hu) (hqF u hu)
  set S := ∑ u ∈ F, q u * Qf u with hS
  have hlhs : ∑ u ∈ F, q u * (Q - Qf u - α * W) = Q * Q - S - α * W * Q := by
    simp only [mul_sub, sum_sub_distrib, ← sum_mul, ← hQ]
    ring
  -- each comparison across folds goes at most one way
  have hpair : 2 * ∑ u ∈ F, q u * a u ≤ Q * Q - S := by
    have hX : ∑ u ∈ F, q u * a u =
        ∑ u ∈ F, ∑ j ∈ F, if f j ≠ f u ∧ s u j < s j u then q u * q j else 0 := by
      refine sum_congr rfl fun u _ => ?_
      rw [ha, mul_sum]
      exact sum_congr rfl fun j _ => by split_ifs <;> ring
    have hY : ∑ u ∈ F, q u * a u =
        ∑ u ∈ F, ∑ j ∈ F, if f u ≠ f j ∧ s j u < s u j then q u * q j else 0 := by
      rw [hX, sum_comm]
      exact sum_congr rfl fun u _ => sum_congr rfl fun j _ => by split_ifs <;> ring
    have hQQ : Q * Q - S = ∑ u ∈ F, ∑ j ∈ F, if f j = f u then 0 else q u * q j := by
      have hall : Q * Q = ∑ u ∈ F, ∑ j ∈ F, q u * q j := by rw [hQ, sum_mul_sum]
      have hS' : S = ∑ u ∈ F, ∑ j ∈ F, if f j = f u then q u * q j else 0 := by
        rw [hS]
        refine sum_congr rfl fun u _ => ?_
        rw [hQf, mul_sum]
        exact sum_congr rfl fun j _ => by split_ifs <;> ring
      rw [hall, hS', ← sum_sub_distrib]
      refine sum_congr rfl fun u _ => ?_
      rw [← sum_sub_distrib]
      exact sum_congr rfl fun j _ => by split_ifs <;> ring
    rw [hQQ, two_mul]
    nth_rw 1 [hX]
    rw [hY, ← sum_add_distrib]
    refine sum_le_sum fun u hu => ?_
    rw [← sum_add_distrib]
    refine sum_le_sum fun j hj => ?_
    have hprod : 0 ≤ q u * q j := mul_nonneg (hqF u hu) (hqF j hj)
    by_cases hfj : f j = f u
    · have hfu : ¬ f u ≠ f j := fun h => h hfj.symm
      simp only [hfj, ne_eq, not_true_eq_false, false_and, if_false, if_true, add_zero, le_refl]
    · have hfu : f u ≠ f j := fun h => hfj h.symm
      rw [if_neg hfj]
      by_cases h1 : s u j < s j u
      · have h2 : ¬ s j u < s u j := fun h2 => lt_asymm h1 h2
        rw [if_pos ⟨hfj, h1⟩, if_neg (fun h => h2 h.2), add_zero]
      · rw [if_neg (fun h => h1 h.2), zero_add]
        split_ifs
        · exact le_rfl
        · exact hprod
  -- the refused units' same-fold weight is at most the heaviest fold's
  have hSle : S ≤ Wmax * Q := by
    rw [hS, hQ, mul_sum]
    refine sum_le_sum fun u hu => ?_
    rw [mul_comm]
    refine mul_le_mul_of_nonneg_right ?_ (hqF u hu)
    calc Qf u = ∑ j ∈ F.filter (fun j => f j = f u), q j := by rw [hQf, sum_filter]
      _ ≤ ∑ j ∈ C.filter (fun j => f j = f u), q j :=
          sum_le_sum_of_subset_of_nonneg (filter_subset_filter _ hFC)
            fun j hj _ => hq j (mem_filter.1 hj).1
      _ ≤ Wmax := hWmax u (hFC hu)
  have hQ0 : 0 ≤ Q := sum_nonneg hqF
  have hW0' : 0 ≤ W := sum_nonneg hq
  have hkey : Q * Q ≤ (2 * α * W + Wmax) * Q := by nlinarith
  rcases hQ0.eq_or_lt with h0 | hpos
  · rw [← h0]
    nlinarith
  · exact le_of_mul_le_mul_right (by linarith) hpos

/-- `ratio_prodPerm` for any real function: averaging over the fold and slot permutations is
averaging over the units. -/
theorem sum_prodPerm_apply {β' γ : Type*} [Fintype β'] [Fintype γ] [DecidableEq β']
    [DecidableEq γ] (t : β' × γ) (g : β' × γ → ℝ) :
    (∑ p : Equiv.Perm β' × Equiv.Perm γ, g (Equiv.prodCongr p.1 p.2 t)) /
        Fintype.card (Equiv.Perm β' × Equiv.Perm γ) =
      (∑ u, g u) / (Fintype.card β' * Fintype.card γ) := by
  set Y : β' → γ → ℝ := fun b c => g (b, c) with hY
  have step1 : ∀ π : Equiv.Perm β', (∑ ρ : Equiv.Perm γ, Y (π t.1) (ρ t.2)) /
      Fintype.card (Equiv.Perm γ) = (∑ c, Y (π t.1) c) / Fintype.card γ :=
    fun π => sum_perm_apply t.2 (fun c => Y (π t.1) c)
  have step2 := sum_perm_apply t.1 (fun b => (∑ c, Y b c) / Fintype.card γ)
  have hcount : ∑ b, ∑ c, Y b c = ∑ u, g u := by rw [← Fintype.sum_prod_type']
  rw [Fintype.card_prod, Fintype.sum_prod_type]
  push_cast
  calc (∑ π : Equiv.Perm β', ∑ ρ : Equiv.Perm γ, g (Equiv.prodCongr (π, ρ).1 (π, ρ).2 t)) /
          (Fintype.card (Equiv.Perm β') * Fintype.card (Equiv.Perm γ))
      = (∑ π : Equiv.Perm β', (∑ ρ : Equiv.Perm γ, Y (π t.1) (ρ t.2)) /
          Fintype.card (Equiv.Perm γ)) / Fintype.card (Equiv.Perm β') := by
        rw [← sum_div, div_div, mul_comm]
        rfl
    _ = (∑ π : Equiv.Perm β', (∑ c, Y (π t.1) c) / Fintype.card γ) /
          Fintype.card (Equiv.Perm β') := by simp only [step1]
    _ = (∑ b, (∑ c, Y b c) / Fintype.card γ) / Fintype.card β' := step2
    _ = (∑ u, g u) / (Fintype.card β' * Fintype.card γ) := by
        rw [← sum_div, div_div, mul_comm, hcount]

/-- **The tournament across folds with likelihood-ratio weights** (our statement, slice 9's CV-b
weighted). Units are `(fold, slot)`; the test unit's fold holds it and further units from the
training law, which the interval never reads, so given the bag the test unit is unit `u` with
chance `w(z_u)/W` and the folds are uniform: the fold and slot permutations, each weighed by the
test unit's `w`. The pair score is equivariant under them (a symmetric fit to the other folds).
Then the weighted cross-fold rule, its total the weight outside the test unit's fold plus its own
(weighted CV+), refuses with at most `2α + W_max/W` of the weight, written here times `W`. -/
theorem coverage_across_likelihoodRatio {β' γ Z α' : Type*} [Fintype β'] [Fintype γ]
    [DecidableEq β'] [DecidableEq γ] [LinearOrder α'] (t : β' × γ)
    (S : (β' × γ → Z) → β' × γ → β' × γ → α')
    (hS : ∀ (p : Equiv.Perm β' × Equiv.Perm γ) z i j, S (z ∘ Equiv.prodCongr p.1 p.2) i j =
      S z (Equiv.prodCongr p.1 p.2 i) (Equiv.prodCongr p.1 p.2 j))
    (z : β' × γ → Z) (w : Z → ℝ) (hw : ∀ v, 0 ≤ w v) {α Wmax : ℝ} (hα0 : 0 ≤ α) (hα1 : α ≤ 1)
    (hW0 : 0 ≤ Wmax)
    (hWmax : ∀ b : β', ∑ j ∈ univ.filter (fun j : β' × γ => j.1 = b), w (z j) ≤ Wmax) :
    (∑ p : Equiv.Perm β' × Equiv.Perm γ,
      (if (1 - α) * (∑ u, w (z (Equiv.prodCongr p.1 p.2 u)) -
          ∑ j ∈ univ.filter (fun j : β' × γ => j.1 = t.1), w (z (Equiv.prodCongr p.1 p.2 j)) +
          w (z (Equiv.prodCongr p.1 p.2 t))) ≤
        strangeAcrossWeight univ Prod.fst (S (z ∘ Equiv.prodCongr p.1 p.2))
          (fun j => w (z (Equiv.prodCongr p.1 p.2 j))) t
        then w (z (Equiv.prodCongr p.1 p.2 t)) else 0)) * ∑ u, w (z u) ≤
      (2 * α * ∑ u, w (z u) + Wmax) *
        ∑ p : Equiv.Perm β' × Equiv.Perm γ, w (z (Equiv.prodCongr p.1 p.2 t)) := by
  classical
  set R : β' × γ → Prop := fun u => (1 - α) * (∑ j, w (z j) -
    ∑ j ∈ univ.filter (fun j : β' × γ => j.1 = u.1), w (z j) + w (z u)) ≤
      strangeAcrossWeight univ Prod.fst (S z) (fun j => w (z j)) u with hR
  set Fn : β' × γ → ℝ := fun u => if R u then w (z u) else 0 with hFn
  -- each permutation: the test unit refused in the permuted data is its image refused in `z`
  have hpt : ∀ p : Equiv.Perm β' × Equiv.Perm γ,
      (if (1 - α) * (∑ u, w (z (Equiv.prodCongr p.1 p.2 u)) -
          ∑ j ∈ univ.filter (fun j : β' × γ => j.1 = t.1), w (z (Equiv.prodCongr p.1 p.2 j)) +
          w (z (Equiv.prodCongr p.1 p.2 t))) ≤
        strangeAcrossWeight univ Prod.fst (S (z ∘ Equiv.prodCongr p.1 p.2))
          (fun j => w (z (Equiv.prodCongr p.1 p.2 j))) t
        then w (z (Equiv.prodCongr p.1 p.2 t)) else 0) = Fn (Equiv.prodCongr p.1 p.2 t) := by
    intro p
    set g := Equiv.prodCongr p.1 p.2 with hg
    have hfst : ∀ i j : β' × γ, (g i).1 = (g j).1 ↔ i.1 = j.1 := fun i j => by
      simp only [hg, Equiv.prodCongr_apply, Prod.map_fst, p.1.injective.eq_iff]
    have htot : ∑ u, w (z (g u)) = ∑ u, w (z u) := Equiv.sum_comp g (fun u => w (z u))
    have hfold : ∑ j ∈ univ.filter (fun j : β' × γ => j.1 = t.1), w (z (g j)) =
        ∑ j ∈ univ.filter (fun j : β' × γ => j.1 = (g t).1), w (z j) :=
      sum_equiv g (fun j => by simp only [mem_filter, mem_univ, true_and, hfst]) fun j _ => rfl
    have hSg : ∀ i j, S (z ∘ g) i j = S z (g i) (g j) := fun i j => hS p z i j
    have hstr : strangeAcrossWeight univ Prod.fst (S (z ∘ g)) (fun j => w (z (g j))) t =
        strangeAcrossWeight univ Prod.fst (S z) (fun j => w (z j)) (g t) := by
      unfold strangeAcrossWeight strangerAcross
      refine sum_equiv g (fun j => ?_) fun j _ => rfl
      simp only [mem_filter, mem_univ, true_and, hSg]
      constructor
      · rintro ⟨h1, h2⟩
        exact ⟨fun h => h1 ((hfst j t).1 h), h2⟩
      · rintro ⟨h1, h2⟩
        exact ⟨fun h => h1 ((hfst j t).2 h), h2⟩
    rw [htot, hfold, hstr]
  rw [Fintype.sum_congr _ _ hpt]
  have hcardG : (0 : ℝ) < Fintype.card (Equiv.Perm β' × Equiv.Perm γ) := by
    exact_mod_cast Fintype.card_pos
  have hcardU : (0 : ℝ) < Fintype.card β' * Fintype.card γ := by
    have : (0 : ℝ) < Fintype.card (β' × γ) := by
      exact_mod_cast Fintype.card_pos_iff.2 ⟨t⟩
    rwa [Fintype.card_prod, Nat.cast_mul] at this
  have hA := sum_prodPerm_apply t Fn
  have hB := sum_prodPerm_apply t (fun u => w (z u))
  rw [div_eq_div_iff hcardG.ne' hcardU.ne'] at hA hB
  -- the counting core on the units of `z`
  have hcount : ∑ u, Fn u ≤ 2 * α * ∑ u, w (z u) + Wmax := by
    have := sum_weight_refused_across_le univ Prod.fst (S z) (fun j => w (z j))
      (fun j _ => hw _) hα0 hα1 hW0 (fun u _ => hWmax u.1)
    rw [sum_filter] at this
    exact this
  have hW : 0 ≤ ∑ u, w (z u) := sum_nonneg fun u _ => hw _
  refine le_of_mul_le_mul_right ?_ hcardU
  calc ((∑ p : Equiv.Perm β' × Equiv.Perm γ, Fn (Equiv.prodCongr p.1 p.2 t)) *
          ∑ u, w (z u)) * (Fintype.card β' * Fintype.card γ)
      = ((∑ p : Equiv.Perm β' × Equiv.Perm γ, Fn (Equiv.prodCongr p.1 p.2 t)) *
          (Fintype.card β' * Fintype.card γ)) * ∑ u, w (z u) := by ring
    _ = ((∑ u, Fn u) * Fintype.card (Equiv.Perm β' × Equiv.Perm γ)) * ∑ u, w (z u) := by
        rw [hA]
    _ ≤ ((2 * α * ∑ u, w (z u) + Wmax) * Fintype.card (Equiv.Perm β' × Equiv.Perm γ)) *
          ∑ u, w (z u) :=
        mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_right hcount hcardG.le) hW
    _ = (2 * α * ∑ u, w (z u) + Wmax) *
          ((∑ u, w (z u)) * Fintype.card (Equiv.Perm β' × Equiv.Perm γ)) := by ring
    _ = (2 * α * ∑ u, w (z u) + Wmax) *
          ((∑ p : Equiv.Perm β' × Equiv.Perm γ, w (z (Equiv.prodCongr p.1 p.2 t))) *
            (Fintype.card β' * Fintype.card γ)) := by rw [hB]
    _ = (2 * α * ∑ u, w (z u) + Wmax) *
          (∑ p : Equiv.Perm β' × Equiv.Perm γ, w (z (Equiv.prodCongr p.1 p.2 t))) *
            (Fintype.card β' * Fintype.card γ) := by ring

end WeightedAcross

/-! ### Slice 14's follow-ups: merged K-fold CV+ under likelihood-ratio weights

With likelihood-ratio weights the folds' weights depend on which unit is the test, so the
fixed-weight average of `WeightedMerge` does not apply. Markov's thresholds `θ W_k`, proportional to
each block's own weight (a fold with the test unit), are fixed given that block as a set, and each
block's weighted rank is super-uniform there: the weighted layer cake. The test unit's weight enters
the refusal only through `(K − 1) w_t`, bounded by any `γ` at least the weight of a unit the rule
can refuse. -/

section WeightedMergeLR

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

omit [Fintype ι] [DecidableEq ι] in
/-- **The weighted layer cake.** Whatever the scores, with weights `q ≥ 0` on `C`, the weight-sum
of `(b − weightLe(u))⁺` is at most `b²/2`: the units whose weight-below is at most `x` weigh at most
`x` (`sum_weight_refused_le`), so `weightLe` is super-uniform under `q`. -/
theorem sum_weight_mul_max_sub_le (C : Finset ι) (s q : ι → ℝ) (hq : ∀ j ∈ C, 0 ≤ q j) (b : ℝ) :
    ∑ u ∈ C, q u * max 0 (b - weightLe C s q u) ≤ b ^ 2 / 2 := by
  have hsu : ∀ x, 0 ≤ x → x < b →
      ∑ u : C, q u * (if weightLe C s q u ≤ x then 1 else 0) ≤ x := by
    intro x hx _
    rw [sum_coe_sort C (fun u => q u * (if weightLe C s q u ≤ x then 1 else 0))]
    simp only [mul_ite, mul_one, mul_zero]
    rw [← sum_filter]
    exact sum_weight_refused_le C s q hq hx
  have h := sum_mul_max_sub_le (fun u : C => q u) (fun u => hq u u.2)
    (fun u : C => weightLe C s q u) (fun u => sum_nonneg fun j hj => hq j (mem_filter.1 hj).1) b
    hsu
  rwa [sum_coe_sort C (fun u => q u * max 0 (b - weightLe C s q u))] at h

/-- **One block's layer cake under likelihood-ratio weights.** Over the reorderings of a block `C`
holding the test position `t` (a fold with the test unit), each weighed by the test unit's `w`, with
a score equivariant on `C` (the fold's fit reads neither): the test unit's `(θ W_C − V)⁺` averages
at most `θ²/2` times `W_C`. Here `W_C` is the block's weight and `V` the block's weight scoring at
most the test unit, itself included. Given the block as a set, the test unit is `u` with chance
`w(z_u)/W_C`, and the threshold `θ W_C` is fixed. -/
theorem sum_permsOn_mul_max_sub_le {Z : Type*} {C : Finset ι} {t : ι} (ht : t ∈ C)
    {S : (ι → Z) → ι → ℝ} (hS : Equivariant C S) (z : ι → Z) (w : Z → ℝ) (hw : ∀ v, 0 ≤ w v)
    (θ : ℝ) :
    ∑ σ ∈ permsOn C, w (z (σ t)) * max 0 (θ * ∑ j ∈ C, w (z (σ j)) -
        weightLe C (S (z ∘ σ)) (fun j => w (z (σ j))) t) ≤
      θ ^ 2 / 2 * ∑ σ ∈ permsOn C, w (z (σ t)) * ∑ j ∈ C, w (z (σ j)) := by
  set WC := ∑ j ∈ C, w (z j) with hWC
  have hsumC : ∀ σ ∈ permsOn C, ∑ j ∈ C, w (z (σ j)) = WC := fun σ hσ =>
    sum_equiv σ (fun j => (mem_iff_of_mem_permsOn hσ j).symm) fun j _ => rfl
  have hwl : ∀ σ ∈ permsOn C, weightLe C (S (z ∘ σ)) (fun j => w (z (σ j))) t =
      weightLe C (S z) (fun j => w (z j)) (σ t) := fun σ hσ => by
    have hfun : S (z ∘ σ) = S z ∘ σ := funext fun i => hS σ hσ z i
    rw [hfun, show (fun j => w (z (σ j))) = (fun j => w (z j)) ∘ σ from rfl, weightLe_comp hσ]
  set f : ι → ℝ := fun u => w (z u) * max 0 (θ * WC - weightLe C (S z) (fun j => w (z j)) u)
    with hf
  have hL : ∑ σ ∈ permsOn C, w (z (σ t)) * max 0 (θ * ∑ j ∈ C, w (z (σ j)) -
      weightLe C (S (z ∘ σ)) (fun j => w (z (σ j))) t) = ∑ σ ∈ permsOn C, f (σ t) :=
    sum_congr rfl fun σ hσ => by rw [hsumC σ hσ, hwl σ hσ, hf]
  have hR : ∑ σ ∈ permsOn C, w (z (σ t)) * ∑ j ∈ C, w (z (σ j)) =
      ∑ σ ∈ permsOn C, (fun u => w (z u) * WC) (σ t) :=
    sum_congr rfl fun σ hσ => by rw [hsumC σ hσ]
  rw [hL, hR]
  have hc : (0 : ℝ) < C.card := by exact_mod_cast card_pos.2 ⟨t, ht⟩
  refine (mul_le_mul_iff_of_pos_left hc).1 ?_
  rw [card_mul_sum_permsOn ht f, mul_left_comm, card_mul_sum_permsOn ht (fun u => w (z u) * WC)]
  have hcake := sum_weight_mul_max_sub_le C (S z) (fun j => w (z j)) (fun j _ => hw _) (θ * WC)
  have hG : (0 : ℝ) ≤ (permsOn C).card := Nat.cast_nonneg _
  calc ((permsOn C).card : ℝ) * ∑ u ∈ C, f u ≤ (permsOn C).card * ((θ * WC) ^ 2 / 2) :=
        mul_le_mul_of_nonneg_left hcake hG
    _ = θ ^ 2 / 2 * ((permsOn C).card * ∑ u ∈ C, w (z u) * WC) := by
        rw [← sum_mul, ← hWC]
        ring

/-- Sums over the blocks: block `C k` is fold `k` with the test position `t`. -/
theorem sum_blocks_filter {κ : Type*} [Fintype κ] {t : ι} (fold : ι → κ)
    {C : κ → Finset ι} (hC : ∀ k i, i ∈ C k ↔ i = t ∨ fold i = k) (f : ι → ℝ)
    (p : κ → ι → Prop) [∀ k, DecidablePred (p k)] :
    ∑ k, ∑ j ∈ (C k).filter (p k), f j =
      ∑ k, (if p k t then f t else 0) +
        ∑ j ∈ (univ.erase t).filter (fun j => p (fold j) j), f j := by
  classical
  have hCk : ∀ k, C k = insert t ((univ.erase t).filter fun j => fold j = k) := fun k => by
    ext j
    by_cases hj : j = t
    · simp [hC, hj]
    · simp [hC, hj]
  have hnot : ∀ k, t ∉ ((univ.erase t).filter fun j => fold j = k).filter (p k) := fun k h => by
    simp at h
  have h1 : ∀ k, ∑ j ∈ (C k).filter (p k), f j = (if p k t then f t else 0) +
      ∑ j ∈ ((univ.erase t).filter fun j => fold j = k).filter (p k), f j := fun k => by
    rw [hCk k, filter_insert]
    split_ifs with h
    · rw [sum_insert (hnot k)]
    · rw [zero_add]
  rw [sum_congr rfl fun k _ => h1 k, sum_add_distrib]
  congr 1
  rw [← sum_fiberwise_of_maps_to (s := (univ.erase t).filter fun j => p (fold j) j) (t := univ)
    (g := fold) fun _ _ => mem_univ _]
  refine sum_congr rfl fun k _ => sum_congr ?_ fun _ _ => rfl
  ext j
  simp only [mem_filter, mem_erase, mem_univ, and_true]
  constructor
  · rintro ⟨⟨hj, hk⟩, hp⟩
    refine ⟨⟨hj, ?_⟩, hk⟩
    rw [hk]
    exact hp
  · rintro ⟨⟨hj, hp⟩, hk⟩
    refine ⟨⟨hj, hk⟩, ?_⟩
    rw [← hk]
    exact hp

omit [DecidableEq ι] in
/-- The blocks' weights: each holds the test position, so they sum to the whole weight plus `K − 1`
copies of the test unit's. -/
theorem sum_blocks {κ : Type*} [Fintype κ] {t : ι} (fold : ι → κ)
    {C : κ → Finset ι} (hC : ∀ k i, i ∈ C k ↔ i = t ∨ fold i = k) (q : ι → ℝ) :
    ∑ k, ∑ j ∈ C k, q j = ∑ i, q i + (Fintype.card κ - 1) * q t := by
  classical
  have h := sum_blocks_filter fold hC q (fun _ _ => True)
  have hf : ∀ k, (C k).filter (fun _ => True) = C k := fun k =>
    filter_true_of_mem fun _ _ => trivial
  have hf' : (univ.erase t).filter (fun _ => True) = univ.erase t :=
    filter_true_of_mem fun _ _ => trivial
  simp only [hf, hf', if_true, sum_const, card_univ, nsmul_eq_mul] at h
  have he := add_sum_erase univ q (mem_univ t)
  rw [h]
  linarith

/-- **A weighted CV+ refusal leaves the blocks' slack large.** Training units in folds `fold`, the
test position `t` in every block `C k`; unit `j`'s comparison with the test unit is under its fold's
score `s (fold j)`. If the units the test unit is stranger than carry at least `1 − α` of the weight
`W`, then for `θ ≤ 1`, `Σ_k (θ W_k − V_k) ≥ (θ − α) W − (1 − θ)(K − 1) γ`. Here `W_k` is block
`k`'s weight, `V_k` its weight scoring at most the test unit, itself included, and `γ` is at least
the test unit's weight when that is at most `α W`, which a refusal needs. -/
theorem blocks_slack_of_refused {κ : Type*} [Fintype κ] {t : ι} (fold : ι → κ)
    {C : κ → Finset ι} (hC : ∀ k i, i ∈ C k ↔ i = t ∨ fold i = k) (s : κ → ι → ℝ) (q : ι → ℝ)
    (hq : ∀ i, 0 ≤ q i) {α θ γ : ℝ} (hθ1 : θ ≤ 1) (hγ : q t ≤ α * ∑ i, q i → q t ≤ γ)
    (href : (1 - α) * ∑ i, q i ≤
      ∑ j ∈ univ.filter (fun j => j ≠ t ∧ s (fold j) t < s (fold j) j), q j) :
    (θ - α) * ∑ i, q i - (1 - θ) * (Fintype.card κ - 1) * γ ≤
      ∑ k, (θ * ∑ j ∈ C k, q j - weightLe (C k) (s k) q t) := by
  have hK : (1 : ℝ) ≤ Fintype.card κ := by exact_mod_cast Fintype.card_pos_iff.2 ⟨fold t⟩
  have hA := sum_blocks fold hC q
  have hB := sum_blocks_filter fold hC q (fun k j => s k j ≤ s k t)
  simp only [le_refl, if_true, sum_const, card_univ, nsmul_eq_mul] at hB
  have hV : ∑ k, weightLe (C k) (s k) q t =
      ∑ k, ∑ j ∈ (C k).filter (fun j => s k j ≤ s k t), q j := rfl
  have he := add_sum_erase univ q (mem_univ t)
  have hsplit := sum_filter_add_sum_filter_not (univ.erase t)
    (fun j => s (fold j) j ≤ s (fold j) t) q
  have hbeat : (univ.erase t).filter (fun j => ¬ s (fold j) j ≤ s (fold j) t) =
      univ.filter (fun j => j ≠ t ∧ s (fold j) t < s (fold j) j) := by
    ext j
    simp only [mem_filter, mem_erase, mem_univ, and_true, true_and, not_le]
  rw [hbeat] at hsplit
  have hU0 : 0 ≤ ∑ j ∈ (univ.erase t).filter (fun j => s (fold j) j ≤ s (fold j) t), q j :=
    sum_nonneg fun j _ => hq j
  have hqt : q t ≤ α * ∑ i, q i := by linarith
  have hqγ := hγ hqt
  have hprod : 0 ≤ (1 - θ) * (Fintype.card κ - 1) * (γ - q t) :=
    mul_nonneg (mul_nonneg (by linarith) (by linarith)) (by linarith)
  rw [sum_sub_distrib, ← mul_sum, hA, hV, hB]
  nlinarith

/-- **Merged K-fold CV+ with likelihood-ratio weights, Markov's form** (our statement). Over the
orderings of a data set `z`, each weighed by the test unit's likelihood ratio `w`, the folds are
given by `fold` and block `C k` is fold `k` with the test position `t`; fold `k`'s score `S k` is
equivariant on its block (a fit to the other folds). Weighted CV+ refuses when the units the test
unit is stranger than, each under its fold's fit, carry at least `1 − α` of the weight `W`, the test
unit's own included. For `θ ≤ 1`, and `γ` at least every weight at most `α W`, the refused weight
times `(θ − α) W − (1 − θ)(K − 1) γ` is at most `θ²/2` times the ordering sum of
`w_t (W + (K − 1) w_t)`. -/
theorem cvPlus_merged_likelihoodRatio_markov {Z κ : Type*} [Fintype κ] {t : ι}
    (fold : ι → κ) {C : κ → Finset ι} (hC : ∀ k i, i ∈ C k ↔ i = t ∨ fold i = k)
    (S : κ → (ι → Z) → ι → ℝ) (hS : ∀ k, Equivariant (C k) (S k)) (z : ι → Z) (w : Z → ℝ)
    (hw : ∀ v, 0 ≤ w v) {α θ γ : ℝ} (hθ1 : θ ≤ 1)
    (hγ : ∀ u, w (z u) ≤ α * ∑ i, w (z i) → w (z u) ≤ γ) :
    ((θ - α) * ∑ i, w (z i) - (1 - θ) * (Fintype.card κ - 1) * γ) *
        ∑ σ : Equiv.Perm ι, (if (1 - α) * ∑ i, w (z (σ i)) ≤ ∑ j ∈ univ.filter (fun j =>
          j ≠ t ∧ S (fold j) (z ∘ σ) t < S (fold j) (z ∘ σ) j), w (z (σ j))
          then w (z (σ t)) else 0) ≤
      θ ^ 2 / 2 * ∑ σ : Equiv.Perm ι, w (z (σ t)) *
        (∑ i, w (z i) + (Fintype.card κ - 1) * w (z (σ t))) := by
  set F : κ → Equiv.Perm ι → ℝ := fun k σ => w (z (σ t)) * max 0 (θ * ∑ j ∈ C k, w (z (σ j)) -
    weightLe (C k) (S k (z ∘ σ)) (fun j => w (z (σ j))) t) with hF
  -- each ordering: a refusal leaves the blocks' slack large
  have hpt : ∀ σ : Equiv.Perm ι, ((θ - α) * ∑ i, w (z i) - (1 - θ) * (Fintype.card κ - 1) * γ) *
      (if (1 - α) * ∑ i, w (z (σ i)) ≤ ∑ j ∈ univ.filter (fun j =>
        j ≠ t ∧ S (fold j) (z ∘ σ) t < S (fold j) (z ∘ σ) j), w (z (σ j))
        then w (z (σ t)) else 0) ≤ ∑ k, F k σ := by
    intro σ
    have hWσ : ∑ i, w (z (σ i)) = ∑ i, w (z i) := Equiv.sum_comp σ (fun i => w (z i))
    have hF0 : 0 ≤ ∑ k, F k σ := sum_nonneg fun k _ => mul_nonneg (hw _) (le_max_left _ _)
    split_ifs with h
    · have hsl := blocks_slack_of_refused (α := α) (θ := θ) (γ := γ)
        (s := fun k => S k (z ∘ σ)) fold hC (fun j => w (z (σ j))) (fun i => hw _) hθ1
        (fun h' => hγ (σ t) (by rw [← hWσ]; exact h')) h
      beta_reduce at hsl
      rw [hWσ] at hsl
      calc ((θ - α) * ∑ i, w (z i) - (1 - θ) * (Fintype.card κ - 1) * γ) * w (z (σ t))
          ≤ (∑ k, (θ * ∑ j ∈ C k, w (z (σ j)) -
              weightLe (C k) (S k (z ∘ σ)) (fun j => w (z (σ j))) t)) * w (z (σ t)) :=
            mul_le_mul_of_nonneg_right hsl (hw _)
        _ ≤ (∑ k, max 0 (θ * ∑ j ∈ C k, w (z (σ j)) -
              weightLe (C k) (S k (z ∘ σ)) (fun j => w (z (σ j))) t)) * w (z (σ t)) :=
            mul_le_mul_of_nonneg_right (sum_le_sum fun k _ => le_max_right _ _) (hw _)
        _ = ∑ k, F k σ := by
            rw [sum_mul]
            refine sum_congr rfl fun k _ => ?_
            simp only [hF]
            ring
    · rw [mul_zero]
      exact hF0
  -- each block: its layer cake over the cosets of its reorderings
  have hre : ∀ (H : Finset (Equiv.Perm ι)) (g : Equiv.Perm ι → ℝ),
      (H.card : ℝ) * ∑ σ : Equiv.Perm ι, g σ = ∑ σ : Equiv.Perm ι, ∑ τ ∈ H, g (σ * τ) := by
    intro H g
    rw [sum_comm, sum_congr rfl fun τ _ =>
      Fintype.sum_equiv (Equiv.mulRight τ) (fun σ => g (σ * τ)) g fun σ => rfl, sum_const,
      nsmul_eq_mul]
  have hblock : ∀ k, ∑ σ : Equiv.Perm ι, F k σ ≤
      θ ^ 2 / 2 * ∑ σ : Equiv.Perm ι, w (z (σ t)) * ∑ j ∈ C k, w (z (σ j)) := by
    intro k
    have htk : t ∈ C k := (hC k t).2 (Or.inl rfl)
    have hH : (0 : ℝ) < (permsOn (C k)).card := by exact_mod_cast card_permsOn_pos (C k)
    refine (mul_le_mul_iff_of_pos_left hH).1 ?_
    rw [hre (permsOn (C k)) (F k), mul_left_comm,
      hre (permsOn (C k)) (fun σ => w (z (σ t)) * ∑ j ∈ C k, w (z (σ j))), mul_sum]
    refine sum_le_sum fun σ _ => ?_
    exact sum_permsOn_mul_max_sub_le htk (hS k) (z ∘ σ) w hw θ
  -- assemble
  have hsum := sum_le_sum fun σ (_ : σ ∈ (univ : Finset (Equiv.Perm ι))) => hpt σ
  rw [← mul_sum] at hsum
  refine hsum.trans ?_
  rw [sum_comm]
  calc ∑ k, ∑ σ : Equiv.Perm ι, F k σ
      ≤ ∑ k, θ ^ 2 / 2 * ∑ σ : Equiv.Perm ι, w (z (σ t)) * ∑ j ∈ C k, w (z (σ j)) :=
        sum_le_sum fun k _ => hblock k
    _ = θ ^ 2 / 2 * ∑ σ : Equiv.Perm ι, w (z (σ t)) * ∑ k, ∑ j ∈ C k, w (z (σ j)) := by
        rw [← mul_sum, sum_comm]
        simp only [mul_sum]
    _ = θ ^ 2 / 2 * ∑ σ : Equiv.Perm ι, w (z (σ t)) *
          (∑ i, w (z i) + (Fintype.card κ - 1) * w (z (σ t))) := by
        congr 1
        refine sum_congr rfl fun σ _ => ?_
        rw [sum_blocks fold hC (fun j => w (z (σ j))), Equiv.sum_comp σ (fun i => w (z i))]

/-- **Merged K-fold CV+ with likelihood-ratio weights** (our statement; follow-up 2's open item).
In the setting of `cvPlus_merged_likelihoodRatio_markov`, with `γ` at least every unit's weight,
weighted CV+ refuses with at most `min(1, 2c(γ))` of the weight, where
`c(γ) = (α W + (K − 1) γ)/(W + (K − 1) γ)`: the position-weight bound `2c` with the test unit's
weight replaced by `γ`. At equal weights it is slice 9's merged bound. -/
theorem cvPlus_merged_likelihoodRatio {Z κ : Type*} [Fintype κ] {t : ι}
    (fold : ι → κ) {C : κ → Finset ι} (hC : ∀ k i, i ∈ C k ↔ i = t ∨ fold i = k)
    (S : κ → (ι → Z) → ι → ℝ) (hS : ∀ k, Equivariant (C k) (S k)) (z : ι → Z) (w : Z → ℝ)
    (hw : ∀ v, 0 ≤ w v) {α γ : ℝ} (hα : 0 < α) (hW : 0 < ∑ i, w (z i))
    (hγ : ∀ u, w (z u) ≤ γ) :
    ∑ σ : Equiv.Perm ι, (if (1 - α) * ∑ i, w (z (σ i)) ≤ ∑ j ∈ univ.filter (fun j =>
        j ≠ t ∧ S (fold j) (z ∘ σ) t < S (fold j) (z ∘ σ) j), w (z (σ j))
        then w (z (σ t)) else 0) ≤
      min 1 (2 * (α * ∑ i, w (z i) + (Fintype.card κ - 1) * γ) /
        (∑ i, w (z i) + (Fintype.card κ - 1) * γ)) * ∑ σ : Equiv.Perm ι, w (z (σ t)) := by
  have hK : (1 : ℝ) ≤ Fintype.card κ := by exact_mod_cast Fintype.card_pos_iff.2 ⟨fold t⟩
  have hγ0 : 0 ≤ γ := (hw _).trans (hγ t)
  have hRT : ∑ σ : Equiv.Perm ι, (if (1 - α) * ∑ i, w (z (σ i)) ≤ ∑ j ∈ univ.filter (fun j =>
      j ≠ t ∧ S (fold j) (z ∘ σ) t < S (fold j) (z ∘ σ) j), w (z (σ j))
      then w (z (σ t)) else 0) ≤ ∑ σ : Equiv.Perm ι, w (z (σ t)) :=
    sum_le_sum fun σ _ => by split_ifs <;> linarith [hw (z (σ t))]
  have hm := cvPlus_merged_likelihoodRatio_markov (α := α) (γ := γ)
    (θ := 2 * (α * ∑ i, w (z i) + (Fintype.card κ - 1) * γ) /
      (∑ i, w (z i) + (Fintype.card κ - 1) * γ)) fold hC S hS z w hw
  set W := ∑ i, w (z i) with hWdef
  set L := α * W + (Fintype.card κ - 1) * γ with hLdef
  set N := W + (Fintype.card κ - 1) * γ with hNdef
  have hN : 0 < N := by nlinarith
  have hL : 0 < L := by nlinarith [mul_pos hα hW]
  rcases le_or_gt N (2 * L) with h2 | h2
  · rw [min_eq_left ((one_le_div hN).2 h2), one_mul]
    exact hRT
  · have hθ1 : 2 * L / N ≤ 1 := (div_le_one hN).2 h2.le
    rw [min_eq_right hθ1]
    have hm' := hm hθ1 (fun u _ => hγ u)
    have hD : (2 * L / N - α) * W - (1 - 2 * L / N) * (Fintype.card κ - 1) * γ = L := by
      have e : (2 * L / N - α) * W - (1 - 2 * L / N) * (Fintype.card κ - 1) * γ =
          2 * L / N * N - L := by
        rw [hLdef, hNdef]
        ring
      rw [e, div_mul_cancel₀ _ hN.ne']
      ring
    rw [hD] at hm'
    have hS2 : ∑ σ : Equiv.Perm ι, w (z (σ t)) * (W + (Fintype.card κ - 1) * w (z (σ t))) ≤
        N * ∑ σ : Equiv.Perm ι, w (z (σ t)) := by
      rw [mul_sum]
      refine sum_le_sum fun σ _ => ?_
      rw [hNdef]
      nlinarith [mul_nonneg (mul_nonneg (hw (z (σ t))) (sub_nonneg.2 hK))
        (sub_nonneg.2 (hγ (σ t)))]
    have hN' : N ≠ 0 := hN.ne'
    refine le_of_mul_le_mul_left ?_ hL
    calc L * ∑ σ : Equiv.Perm ι, (if (1 - α) * ∑ i, w (z (σ i)) ≤ ∑ j ∈ univ.filter (fun j =>
          j ≠ t ∧ S (fold j) (z ∘ σ) t < S (fold j) (z ∘ σ) j), w (z (σ j))
          then w (z (σ t)) else 0)
        ≤ (2 * L / N) ^ 2 / 2 * ∑ σ : Equiv.Perm ι, w (z (σ t)) *
            (W + (Fintype.card κ - 1) * w (z (σ t))) := hm'
      _ ≤ (2 * L / N) ^ 2 / 2 * (N * ∑ σ : Equiv.Perm ι, w (z (σ t))) :=
          mul_le_mul_of_nonneg_left hS2 (by positivity)
      _ = L * (2 * L / N * ∑ σ : Equiv.Perm ι, w (z (σ t))) := by
          field_simp

/-- **The refined form**: with `γ` at least every weight at most `α W` (the units the rule can
refuse) and `2L ≤ N`, `L = α W + (K − 1) γ`, `N = W + (K − 1) γ`, the refused weight is at most
`2L (W + (K − 1) e)/N²` of the whole, `e = Σ w²/W` the test unit's mean weight. Written here times
`N² W`. It reads `γ` below the bag's largest weight when that unit is too heavy to refuse. -/
theorem cvPlus_merged_likelihoodRatio_refined {Z κ : Type*} [Fintype κ] {t : ι}
    (fold : ι → κ) {C : κ → Finset ι} (hC : ∀ k i, i ∈ C k ↔ i = t ∨ fold i = k)
    (S : κ → (ι → Z) → ι → ℝ) (hS : ∀ k, Equivariant (C k) (S k)) (z : ι → Z) (w : Z → ℝ)
    (hw : ∀ v, 0 ≤ w v) {α γ : ℝ} (hα : 0 < α) (hγ0 : 0 ≤ γ) (hW : 0 < ∑ i, w (z i))
    (hγ : ∀ u, w (z u) ≤ α * ∑ i, w (z i) → w (z u) ≤ γ)
    (h2 : 2 * (α * ∑ i, w (z i) + (Fintype.card κ - 1) * γ) ≤
      ∑ i, w (z i) + (Fintype.card κ - 1) * γ) :
    (∑ σ : Equiv.Perm ι, (if (1 - α) * ∑ i, w (z (σ i)) ≤ ∑ j ∈ univ.filter (fun j =>
        j ≠ t ∧ S (fold j) (z ∘ σ) t < S (fold j) (z ∘ σ) j), w (z (σ j))
        then w (z (σ t)) else 0)) * (∑ i, w (z i) + (Fintype.card κ - 1) * γ) ^ 2 *
        ∑ i, w (z i) ≤
      2 * (α * ∑ i, w (z i) + (Fintype.card κ - 1) * γ) *
        ((∑ i, w (z i)) ^ 2 + (Fintype.card κ - 1) * ∑ i, w (z i) ^ 2) *
          ∑ σ : Equiv.Perm ι, w (z (σ t)) := by
  have hK : (1 : ℝ) ≤ Fintype.card κ := by exact_mod_cast Fintype.card_pos_iff.2 ⟨fold t⟩
  have hm := cvPlus_merged_likelihoodRatio_markov (α := α) (γ := γ)
    (θ := 2 * (α * ∑ i, w (z i) + (Fintype.card κ - 1) * γ) /
      (∑ i, w (z i) + (Fintype.card κ - 1) * γ)) fold hC S hS z w hw
  set W := ∑ i, w (z i) with hWdef
  set L := α * W + (Fintype.card κ - 1) * γ with hLdef
  set N := W + (Fintype.card κ - 1) * γ with hNdef
  set R := ∑ σ : Equiv.Perm ι, (if (1 - α) * ∑ i, w (z (σ i)) ≤ ∑ j ∈ univ.filter (fun j =>
      j ≠ t ∧ S (fold j) (z ∘ σ) t < S (fold j) (z ∘ σ) j), w (z (σ j))
      then w (z (σ t)) else 0) with hRdef
  set X := ∑ σ : Equiv.Perm ι, w (z (σ t)) * (W + (Fintype.card κ - 1) * w (z (σ t))) with hXdef
  set T := ∑ σ : Equiv.Perm ι, w (z (σ t)) with hTdef
  have hN : 0 < N := by nlinarith
  have hL : 0 < L := by nlinarith [mul_pos hα hW]
  have hθ1 : 2 * L / N ≤ 1 := (div_le_one hN).2 h2
  have hm' := hm hθ1 hγ
  have hD : (2 * L / N - α) * W - (1 - 2 * L / N) * (Fintype.card κ - 1) * γ = L := by
    have e : (2 * L / N - α) * W - (1 - 2 * L / N) * (Fintype.card κ - 1) * γ =
        2 * L / N * N - L := by
      rw [hLdef, hNdef]
      ring
    rw [e, div_mul_cancel₀ _ hN.ne']
    ring
  rw [hD] at hm'
  have hN' : N ≠ 0 := hN.ne'
  have hRN : R * N ^ 2 ≤ 2 * L * X := by
    have key : L * (R * N ^ 2) ≤ L * (2 * L * X) := by
      calc L * (R * N ^ 2) = (L * R) * N ^ 2 := by ring
        _ ≤ ((2 * L / N) ^ 2 / 2 * X) * N ^ 2 := mul_le_mul_of_nonneg_right hm' (sq_nonneg N)
        _ = L * (2 * L * X) := by
            field_simp
    exact le_of_mul_le_mul_left key hL
  -- the ordering averages: given the bag, the test unit is `u` with chance `w(z_u)/W`
  have hG : (0 : ℝ) < Fintype.card (Equiv.Perm ι) := by exact_mod_cast Fintype.card_pos
  have hI : (0 : ℝ) < Fintype.card ι := by exact_mod_cast Fintype.card_pos_iff.2 ⟨t⟩
  have hA := sum_perm_apply t (fun u => w (z u) * (W + (Fintype.card κ - 1) * w (z u)))
  have hB := sum_perm_apply t (fun u => w (z u))
  rw [div_eq_div_iff hG.ne' hI.ne'] at hA hB
  rw [← hXdef] at hA
  rw [← hTdef, ← hWdef] at hB
  have hsq : ∑ u, w (z u) * (W + (Fintype.card κ - 1) * w (z u)) =
      W ^ 2 + (Fintype.card κ - 1) * ∑ i, w (z i) ^ 2 := by
    have e : ∀ u, w (z u) * (W + (Fintype.card κ - 1) * w (z u)) =
        w (z u) * W + (Fintype.card κ - 1) * w (z u) ^ 2 := fun u => by ring
    rw [sum_congr rfl fun u _ => e u, sum_add_distrib, ← sum_mul, ← mul_sum, ← hWdef]
    ring
  rw [hsq] at hA
  have hXW : X * W = (W ^ 2 + (Fintype.card κ - 1) * ∑ i, w (z i) ^ 2) * T := by
    have h3 : X * W * Fintype.card ι =
        (W ^ 2 + (Fintype.card κ - 1) * ∑ i, w (z i) ^ 2) * T * Fintype.card ι := by
      calc X * W * Fintype.card ι = (X * Fintype.card ι) * W := by ring
        _ = ((W ^ 2 + (Fintype.card κ - 1) * ∑ i, w (z i) ^ 2) *
              Fintype.card (Equiv.Perm ι)) * W := by rw [hA]
        _ = (W ^ 2 + (Fintype.card κ - 1) * ∑ i, w (z i) ^ 2) *
              (W * Fintype.card (Equiv.Perm ι)) := by ring
        _ = (W ^ 2 + (Fintype.card κ - 1) * ∑ i, w (z i) ^ 2) * (T * Fintype.card ι) := by
            rw [hB]
        _ = (W ^ 2 + (Fintype.card κ - 1) * ∑ i, w (z i) ^ 2) * T * Fintype.card ι := by ring
    exact mul_right_cancel₀ hI.ne' h3
  calc R * N ^ 2 * W ≤ 2 * L * X * W := mul_le_mul_of_nonneg_right hRN hW.le
    _ = 2 * L * (W ^ 2 + (Fintype.card κ - 1) * ∑ i, w (z i) ^ 2) * T := by
        rw [mul_assoc, hXW]
        ring

end WeightedMergeLR

/-! ### Slice 15: posterior conformal prediction (Zhang and Candès, arXiv 2409.19712) -/

section PosteriorConformal

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- **Posterior conformal prediction, given the draw** (their Theorem 1), per data set. The test
unit's value `v` sends a draw `L` with chance `κ v L` (their multinomial `L*` of `π(X_{n+1})`).
Given the data set and `L`, the test unit is unit `u` with chance proportional to `κ (z_u) L`,
their posterior weights `w*`, and the weighted rule with those weights refuses it with at most `α`
of that chance. -/
theorem coverage_pcp {Z Λ : Type*} {C : Finset ι} {t : ι} (ht : t ∈ C) {S : (ι → Z) → ι → ℝ}
    (hS : Equivariant C S) (z : ι → Z) (κ : Z → Λ → ℝ) (hκ : ∀ v L, 0 ≤ κ v L) {α : ℝ}
    (hα : 0 ≤ α) (L : Λ) :
    ∑ σ ∈ permsOn C, (if α < weightLe C (S (z ∘ σ))
        (fun j => κ (z (σ j)) L / ∑ i ∈ C, κ (z i) L) t then 0 else κ (z (σ t)) L) ≤
      α * ∑ σ ∈ permsOn C, κ (z (σ t)) L :=
  coverage_weighted ht hS z (fun v => κ v L) (fun v => hκ v L) hα

/-- **And marginally**: summed over the draws, the refused chance is at most `α`. -/
theorem coverage_pcp_marginal {Z Λ : Type*} [Fintype Λ] {C : Finset ι} {t : ι} (ht : t ∈ C)
    {S : (ι → Z) → ι → ℝ} (hS : Equivariant C S) (z : ι → Z) (κ : Z → Λ → ℝ)
    (hκ : ∀ v L, 0 ≤ κ v L) (hκ1 : ∀ v, ∑ L, κ v L = 1) {α : ℝ} (hα : 0 ≤ α) :
    ∑ L, ∑ σ ∈ permsOn C, (if α < weightLe C (S (z ∘ σ))
        (fun j => κ (z (σ j)) L / ∑ i ∈ C, κ (z i) L) t then 0 else κ (z (σ t)) L) ≤
      α * (permsOn C).card := by
  calc ∑ L, ∑ σ ∈ permsOn C, (if α < weightLe C (S (z ∘ σ))
        (fun j => κ (z (σ j)) L / ∑ i ∈ C, κ (z i) L) t then 0 else κ (z (σ t)) L)
      ≤ ∑ L, α * ∑ σ ∈ permsOn C, κ (z (σ t)) L :=
        sum_le_sum fun L _ => coverage_pcp ht hS z κ hκ hα L
    _ = α * (permsOn C).card := by
        rw [← mul_sum, sum_comm]
        simp only [hκ1, sum_const, nsmul_eq_mul, mul_one]

/-- The law of independent units, unit `i` with law `m i`. -/
def prodLaw {V : Type*} (m : ι → V → ℝ) (z : ι → V) : ℝ := ∏ i, m i (z i)

/-- Independent units: summing a product of per-unit factors over every data set. -/
theorem sum_prodLaw {V : Type*} [Fintype V] (c : ι → V → ℝ) :
    ∑ z : ι → V, ∏ i, c i (z i) = ∏ i, ∑ v, c i v := by
  have h := Finset.prod_univ_sum (t := fun _ : ι => (univ : Finset V)) (f := c)
  rw [Fintype.piFinset_univ] at h
  exact h.symm

/-- Splitting a product over every unit at two distinct units. -/
theorem prod_split_two (x : ι → ℝ) {K t : ι} (hKt : K ≠ t) :
    ∏ i, x i = x K * (x t * ∏ i ∈ (univ.erase K).erase t, x i) := by
  rw [← mul_prod_erase univ _ (mem_univ K),
    ← mul_prod_erase (univ.erase K) _ (mem_erase.2 ⟨Ne.symm hKt, mem_univ t⟩)]

omit [Fintype ι] [DecidableEq ι] in
theorem max_sub_add_max_sub (a b : ℝ) : max 0 (b - a) + max 0 (a - b) = |a - b| := by
  rcases le_total a b with h | h
  · rw [max_eq_right (sub_nonneg.2 h), max_eq_left (sub_nonpos.2 h),
      abs_of_nonpos (sub_nonpos.2 h)]
    ring
  · rw [max_eq_left (sub_nonpos.2 h), max_eq_right (sub_nonneg.2 h),
      abs_of_nonneg (sub_nonneg.2 h)]
    ring

/-- **Swapping two independent units costs at most the L1 distance of their laws**, `2 TV`
(Barber et al.'s Lemma 1, without its square term). -/
theorem swap_tv_prodLaw_le {V : Type*} [Fintype V] (m : ι → V → ℝ) (hm0 : ∀ i v, 0 ≤ m i v)
    (hm1 : ∀ i, ∑ v, m i v = 1) (K t : ι) :
    ∑ z, max 0 (prodLaw m (z ∘ Equiv.swap K t) - prodLaw m z) ≤ ∑ v, |m K v - m t v| := by
  by_cases hKt : K = t
  · subst hKt
    simp only [Equiv.swap_self, Equiv.coe_refl, Function.comp_id, sub_self, max_self,
      sum_const_zero]
    exact sum_nonneg fun v _ => abs_nonneg _
  have hR : ∀ z : ι → V, 0 ≤ ∏ i ∈ (univ.erase K).erase t, m i (z i) :=
    fun z => prod_nonneg fun i _ => hm0 i (z i)
  have hrest : ∀ (c : ι → V → ℝ) (z : ι → V),
      (∀ i, i ≠ K → i ≠ t → c i = m i) →
      ∏ i ∈ (univ.erase K).erase t, c i (z i) = ∏ i ∈ (univ.erase K).erase t, m i (z i) :=
    fun c z hc => prod_congr rfl fun i hi => by
      rw [mem_erase, mem_erase] at hi
      rw [hc i hi.2.1 hi.1]
  have rest_swap : ∀ z : ι → V, ∏ i ∈ (univ.erase K).erase t, m i ((z ∘ Equiv.swap K t) i) =
      ∏ i ∈ (univ.erase K).erase t, m i (z i) := fun z => prod_congr rfl fun i hi => by
    rw [mem_erase, mem_erase] at hi
    simp [Function.comp_apply, Equiv.swap_apply_of_ne_of_ne hi.2.1 hi.1]
  let g₁ : ι → V → ℝ := fun i v =>
    if i = K then max 0 (m t v - m K v) else if i = t then m K v else m i v
  let g₂ : ι → V → ℝ := fun i v =>
    if i = K then m K v else if i = t then max 0 (m K v - m t v) else m i v
  have hg₁ : ∀ i, i ≠ K → i ≠ t → g₁ i = m i := fun i h1 h2 => by
    funext v
    simp [g₁, h1, h2]
  have hg₂ : ∀ i, i ≠ K → i ≠ t → g₂ i = m i := fun i h1 h2 => by
    funext v
    simp [g₂, h1, h2]
  have hpoint : ∀ z : ι → V, max 0 (prodLaw m (z ∘ Equiv.swap K t) - prodLaw m z) ≤
      ∏ i, g₁ i (z i) + ∏ i, g₂ i (z i) := by
    intro z
    rw [prodLaw, prodLaw, prod_split_two (fun i => m i ((z ∘ Equiv.swap K t) i)) hKt,
      prod_split_two (fun i => m i (z i)) hKt, prod_split_two (fun i => g₁ i (z i)) hKt,
      prod_split_two (fun i => g₂ i (z i)) hKt, rest_swap, hrest g₁ z hg₁, hrest g₂ z hg₂]
    simp only [Function.comp_apply, Equiv.swap_apply_left, Equiv.swap_apply_right, g₁, g₂,
      if_neg (Ne.symm hKt), ↓reduceIte]
    have hA := hm0 K (z t)
    have hB := hm0 K (z K)
    refine max_le (add_nonneg (mul_nonneg (le_max_left _ _) (mul_nonneg hA (hR z)))
      (mul_nonneg hB (mul_nonneg (le_max_left _ _) (hR z)))) ?_
    have q1 := mul_le_mul_of_nonneg_left (le_max_right 0 (m t (z K) - m K (z K))) hA
    have q2 := mul_le_mul_of_nonneg_left (le_max_right 0 (m K (z t) - m t (z t))) hB
    nlinarith [mul_le_mul_of_nonneg_right q1 (hR z), mul_le_mul_of_nonneg_right q2 (hR z)]
  have hone : ∀ (c : ι → V → ℝ), (∀ i, i ≠ K → i ≠ t → c i = m i) →
      ∏ i ∈ (univ.erase K).erase t, ∑ v, c i v = 1 := fun c hc =>
    prod_eq_one fun i hi => by
      rw [mem_erase, mem_erase] at hi
      rw [hc i hi.2.1 hi.1, hm1]
  have s1 : ∑ z : ι → V, ∏ i, g₁ i (z i) = ∑ v, max 0 (m t v - m K v) := by
    rw [sum_prodLaw, prod_split_two (fun i => ∑ v, g₁ i v) hKt, hone g₁ hg₁]
    simp only [g₁, if_neg (Ne.symm hKt), ↓reduceIte, hm1, mul_one]
  have s2 : ∑ z : ι → V, ∏ i, g₂ i (z i) = ∑ v, max 0 (m K v - m t v) := by
    rw [sum_prodLaw, prod_split_two (fun i => ∑ v, g₂ i v) hKt, hone g₂ hg₂]
    simp only [g₂, if_neg (Ne.symm hKt), ↓reduceIte, hm1, one_mul, mul_one]
  calc ∑ z, max 0 (prodLaw m (z ∘ Equiv.swap K t) - prodLaw m z)
      ≤ ∑ z : ι → V, (∏ i, g₁ i (z i) + ∏ i, g₂ i (z i)) := sum_le_sum fun z _ => hpoint z
    _ = ∑ v, |m K v - m t v| := by
        rw [sum_add_distrib, s1, s2, ← sum_add_distrib]
        exact sum_congr rfl fun v _ => max_sub_add_max_sub (m K v) (m t v)

omit [Fintype ι] [DecidableEq ι] in
/-- Two mixtures of the same components are at most the L1 distance of their weights apart. -/
theorem mixture_l1_le {V Kc : Type*} [Fintype V] [Fintype Kc] (f : Kc → V → ℝ)
    (hf0 : ∀ k v, 0 ≤ f k v) (hf1 : ∀ k, ∑ v, f k v = 1) (a b : Kc → ℝ) :
    ∑ v, |∑ k, a k * f k v - ∑ k, b k * f k v| ≤ ∑ k, |a k - b k| := by
  calc ∑ v, |∑ k, a k * f k v - ∑ k, b k * f k v|
      = ∑ v, |∑ k, (a k - b k) * f k v| := by
        refine sum_congr rfl fun v _ => ?_
        rw [← sum_sub_distrib]
        exact congrArg _ (sum_congr rfl fun k _ => by ring)
    _ ≤ ∑ v, ∑ k, |a k - b k| * f k v := sum_le_sum fun v _ =>
        (abs_sum_le_sum_abs _ _).trans (le_of_eq (sum_congr rfl fun k _ => by
          rw [abs_mul, abs_of_nonneg (hf0 k v)]))
    _ = ∑ k, |a k - b k| := by
        rw [sum_comm]
        exact sum_congr rfl fun k _ => by rw [← mul_sum, hf1, mul_one]

/-- **PCP's per-cluster gap** (Zhang and Candès, Theorem 3), on a finite space, with half their
constant. Given the covariates the residuals are independent, unit `i`'s law the mixture
`Σ_k π_k(i) f_k`; the weights are fixed, nonnegative, summing to one, the test's the largest; the
score is symmetric. Then the miscoverage is at most `α + Σ_K q_K Σ_k |π_k(K) − π_k(t)|`. They
state `α + 2 Σ_K q_K Σ_k |π_k(K) − π_k(t)|`; through Barber et al.'s swap bound the swap costs
the L1 distance of the two units' laws, at most that of their weights, once. -/
theorem pcp_cluster_gap {V Kc : Type*} [Fintype V] [Fintype Kc] {S : (ι → V) → ι → ℝ}
    (hS : Equivariant univ S) (f : Kc → V → ℝ) (hf0 : ∀ k v, 0 ≤ f k v)
    (hf1 : ∀ k, ∑ v, f k v = 1) (π : ι → Kc → ℝ) (hπ0 : ∀ i k, 0 ≤ π i k)
    (hπ1 : ∀ i, ∑ k, π i k = 1) (q : ι → ℝ) (hq0 : ∀ i, 0 ≤ q i) (hq1 : ∑ i, q i = 1) (t : ι)
    (hqt : ∀ i, q i ≤ q t) {α : ℝ} (hα : 0 ≤ α) :
    ∑ z, prodLaw (fun i v => ∑ k, π i k * f k v) z *
        (if weightLe univ (S z) q t ≤ α then 1 else 0) ≤
      α + ∑ K, q K * ∑ k, |π K k - π t k| := by
  have hm0 : ∀ i v, 0 ≤ ∑ k, π i k * f k v :=
    fun i v => sum_nonneg fun k _ => mul_nonneg (hπ0 i k) (hf0 k v)
  have hm1 : ∀ i, ∑ v, ∑ k, π i k * f k v = 1 := fun i => by
    rw [sum_comm]
    simp only [← mul_sum, hf1, mul_one, hπ1]
  have hP1 : ∑ z, prodLaw (fun i v => ∑ k, π i k * f k v) z = 1 := by
    simp only [prodLaw]
    rw [sum_prodLaw (fun i v => ∑ k, π i k * f k v)]
    simp only [hm1, prod_const_one]
  have h := refused_le_swap_tv hS (prodLaw fun i v => ∑ k, π i k * f k v)
    (fun z => prod_nonneg fun i _ => hm0 i (z i)) q hq0 hq1 t hqt hα
  rw [hP1, mul_one] at h
  refine h.trans (add_le_add le_rfl (sum_le_sum fun K _ =>
    mul_le_mul_of_nonneg_left ?_ (hq0 K)))
  exact (swap_tv_prodLaw_le _ hm0 hm1 K t).trans (mixture_l1_le f hf0 hf1 (π K) (π t))

end PosteriorConformal

/-! ### Slice 16: batch predictive inference

Lee, Tchetgen Tchetgen and Dobriban, arXiv 2409.13990. -/

section Batch

variable {n m : ℕ}

/-- The calibration rank of the `j`-th smallest test unit, 0-based: the calibration units (outside
the test set `A`) below it. Units are their score ranks among the `n + m`, scores distinct. -/
noncomputable def testRank (A : Finset (Fin (n + m))) (j : Fin m) : ℕ :=
  if h : A.card = m then ((univ \ A).filter fun p => p < A.orderEmbOfFin h j).card else 0

/-- The rank sequences: nondecreasing, each at most `n` (their `H`, shifted to 0-based). -/
noncomputable def rankSeqs (n m : ℕ) : Finset (Fin m → ℕ) :=
  (Fintype.piFinset fun _ : Fin m => range (n + 1)).filter Monotone

/-- The units of a set below its `j`-th smallest are its `j` smallest. -/
theorem card_filter_lt_orderEmbOfFin {N : ℕ} (A : Finset (Fin N)) {k : ℕ} (h : A.card = k)
    (j : Fin k) : (A.filter fun p => p < A.orderEmbOfFin h j).card = j := by
  have e : (A.filter fun p => p < A.orderEmbOfFin h j) =
      (Finset.Iio j).map (A.orderEmbOfFin h).toEmbedding := by
    ext p
    simp only [mem_filter, mem_map, Finset.mem_Iio, RelEmbedding.coe_toEmbedding]
    constructor
    · rintro ⟨hp, hlt⟩
      have : p ∈ Set.range (A.orderEmbOfFin h) := by
        rw [Finset.range_orderEmbOfFin]
        exact hp
      obtain ⟨k', rfl⟩ := this
      exact ⟨k', (A.orderEmbOfFin h).lt_iff_lt.1 hlt, rfl⟩
    · rintro ⟨k', hk', rfl⟩
      exact ⟨A.orderEmbOfFin_mem h k', (A.orderEmbOfFin h).lt_iff_lt.2 hk'⟩
  rw [e, card_map, Fin.card_Iio]

theorem card_univ_filter_lt {N : ℕ} (q : Fin N) :
    ((univ : Finset (Fin N)).filter fun p => p < q).card = q := by
  rw [show ((univ : Finset (Fin N)).filter fun p => p < q) = Finset.Iio q by ext; simp]
  exact Fin.card_Iio _

/-- The calibration rank is the test unit's full rank less its rank among the tests. -/
theorem testRank_eq (A : Finset (Fin (n + m))) (h : A.card = m) (j : Fin m) :
    testRank A j = (A.orderEmbOfFin h j : ℕ) - j := by
  rw [testRank, dif_pos h]
  have h1 := card_filter_lt_orderEmbOfFin A h j
  have h2 := card_univ_filter_lt (A.orderEmbOfFin h j)
  have hdisj : Disjoint (A.filter fun p => p < A.orderEmbOfFin h j)
      ((univ \ A).filter fun p => p < A.orderEmbOfFin h j) := by
    rw [disjoint_left]
    intro p hp hq
    simp only [mem_filter, mem_sdiff] at hp hq
    exact hq.1.2 hp.1
  have hunion : (univ : Finset (Fin (n + m))).filter (fun p => p < A.orderEmbOfFin h j) =
      (A.filter fun p => p < A.orderEmbOfFin h j) ∪
        ((univ \ A).filter fun p => p < A.orderEmbOfFin h j) := by
    ext p
    simp only [mem_filter, mem_union, mem_sdiff, mem_univ, true_and]
    tauto
  have := card_union_of_disjoint hdisj
  rw [← hunion, h2, h1] at this
  omega

theorem le_orderEmbOfFin (A : Finset (Fin (n + m))) (h : A.card = m) (j : Fin m) :
    (j : ℕ) ≤ A.orderEmbOfFin h j := by
  calc (j : ℕ) = (A.filter fun p => p < A.orderEmbOfFin h j).card :=
        (card_filter_lt_orderEmbOfFin A h j).symm
    _ ≤ ((univ : Finset (Fin (n + m))).filter fun p => p < A.orderEmbOfFin h j).card :=
        card_le_card fun p hp => by
          simp only [mem_filter] at hp ⊢
          exact ⟨mem_univ _, hp.2⟩
    _ = A.orderEmbOfFin h j := card_univ_filter_lt _

/-- The test set from a rank sequence: the units `r_j + j`. -/
def ofRanks (r : Fin m → ℕ) : Finset (Fin (n + m)) :=
  univ.filter fun p => ∃ j : Fin m, (p : ℕ) = r j + j

theorem strictMono_add_index {r : Fin m → ℕ} (hr : Monotone r) :
    StrictMono fun j : Fin m => r j + (j : ℕ) := fun i j hij => by
  change r i + (i : ℕ) < r j + j
  have := hr hij.le
  have : (i : ℕ) < j := hij
  omega

theorem ofRanks_eq_image {r : Fin m → ℕ} (hr : ∀ j, r j ≤ n) :
    (ofRanks r : Finset (Fin (n + m))) = univ.image fun j : Fin m =>
      (⟨r j + j, by have := hr j; have := j.isLt; omega⟩ : Fin (n + m)) := by
  ext p
  simp only [ofRanks, mem_filter, mem_univ, true_and, mem_image]
  constructor
  · rintro ⟨j, hj⟩
    exact ⟨j, Fin.ext hj.symm⟩
  · rintro ⟨j, rfl⟩
    exact ⟨j, rfl⟩

theorem card_ofRanks {r : Fin m → ℕ} (hr : ∀ j, r j ≤ n) (hm : Monotone r) :
    (ofRanks r : Finset (Fin (n + m))).card = m := by
  rw [ofRanks_eq_image hr, card_image_of_injective _ fun i j hij =>
    (strictMono_add_index hm).injective (congrArg Fin.val hij), card_univ, Fintype.card_fin]

theorem testRank_ofRanks {r : Fin m → ℕ} (hr : ∀ j, r j ≤ n) (hm : Monotone r) :
    testRank (ofRanks r : Finset (Fin (n + m))) = r := by
  have hc := card_ofRanks (n := n) hr hm
  funext j
  rw [testRank_eq _ hc]
  have hu := Finset.orderEmbOfFin_unique hc
    (f := fun j : Fin m => (⟨r j + j, by have := hr j; have := j.isLt; omega⟩ : Fin (n + m)))
    (fun j => by rw [ofRanks_eq_image hr]; exact mem_image_of_mem _ (mem_univ j))
    (fun i j hij => strictMono_add_index hm hij)
  rw [← congrFun hu j]
  simp

theorem ofRanks_testRank (A : Finset (Fin (n + m))) (h : A.card = m) :
    ofRanks (testRank A) = A := by
  ext p
  simp only [ofRanks, mem_filter, mem_univ, true_and]
  constructor
  · rintro ⟨j, hj⟩
    rw [testRank_eq A h, Nat.sub_add_cancel (le_orderEmbOfFin A h j)] at hj
    rw [show p = A.orderEmbOfFin h j from Fin.ext hj]
    exact A.orderEmbOfFin_mem h j
  · intro hp
    have : p ∈ Set.range (A.orderEmbOfFin h) := by
      rw [Finset.range_orderEmbOfFin]
      exact hp
    obtain ⟨j, rfl⟩ := this
    exact ⟨j, by rw [testRank_eq A h, Nat.sub_add_cancel (le_orderEmbOfFin A h j)]⟩

/-- The calibration ranks of a test set form a rank sequence. -/
theorem testRank_mem_rankSeqs (A : Finset (Fin (n + m))) (h : A.card = m) :
    testRank A ∈ rankSeqs n m := by
  rw [rankSeqs, mem_filter, Fintype.mem_piFinset]
  refine ⟨fun j => ?_, fun i j hij => ?_⟩
  · rw [mem_range]
    simp only [testRank, dif_pos h]
    have hcA : (univ \ A).card = n := by
      rw [card_sdiff_of_subset (subset_univ A), card_univ, Fintype.card_fin, h]
      omega
    have := card_le_card (filter_subset (fun p => p < A.orderEmbOfFin h j) (univ \ A))
    rw [hcA] at this
    omega
  · simp only [testRank, dif_pos h]
    exact card_le_card fun p hp => by
      simp only [mem_filter] at hp ⊢
      exact ⟨hp.1, lt_of_lt_of_le hp.2 ((A.orderEmbOfFin h).monotone hij)⟩

/-- **The ranks are a bijection from test sets onto the rank sequences** (their map `M`, with
inverse `r ↦ r_j + j`). With the test set uniform over the `C(n + m, m)` subsets the rank sequence
is uniform over the sequences, `|H| = C(n + m, m)`: the test sets ranked in any `I` number
`|H ∩ I|`. -/
theorem card_testRank_mem (I : Finset (Fin m → ℕ)) :
    ((powersetCard m (univ : Finset (Fin (n + m)))).filter fun A => testRank A ∈ I).card =
      (rankSeqs n m ∩ I).card := by
  refine card_nbij' (fun A => testRank A) (fun r => ofRanks r) ?_ ?_ ?_ ?_
  · intro A hA
    rw [mem_coe, mem_filter, mem_powersetCard] at hA
    exact mem_inter.2 ⟨testRank_mem_rankSeqs A hA.1.2, hA.2⟩
  · intro r hr
    rw [mem_coe, mem_inter, rankSeqs, mem_filter, Fintype.mem_piFinset] at hr
    have hb : ∀ j, r j ≤ n := fun j => Nat.lt_succ_iff.1 (mem_range.1 (hr.1.1 j))
    rw [mem_coe, mem_filter, mem_powersetCard, testRank_ofRanks hb hr.1.2]
    exact ⟨⟨subset_univ _, card_ofRanks hb hr.1.2⟩, hr.2⟩
  · intro A hA
    rw [mem_coe, mem_filter, mem_powersetCard] at hA
    exact ofRanks_testRank A hA.1.2
  · intro r hr
    rw [mem_coe, mem_inter, rankSeqs, mem_filter, Fintype.mem_piFinset] at hr
    exact testRank_ofRanks (fun j => Nat.lt_succ_iff.1 (mem_range.1 (hr.1.1 j))) hr.1.2

/-- `|H| = C(n + m, m)`: the rank sequences number the test sets. -/
theorem card_rankSeqs : (rankSeqs n m).card = (n + m).choose m := by
  have := card_testRank_mem (n := n) (m := m) (rankSeqs n m)
  rw [inter_self, filter_true_of_mem fun A hA =>
    testRank_mem_rankSeqs A (mem_powersetCard.1 hA).2, card_powersetCard, card_univ,
    Fintype.card_fin] at this
  exact this.symm

/-- **The order-statistic sandwich, from above**: the `j`-th smallest test unit lies below the
calibration unit at its calibration rank, when there is one. -/
theorem orderEmbOfFin_lt_cal (A : Finset (Fin (n + m))) (h : A.card = m)
    (hc : (univ \ A).card = n) (j : Fin m) (hlt : testRank A j < n) :
    A.orderEmbOfFin h j < (univ \ A).orderEmbOfFin hc ⟨testRank A j, hlt⟩ := by
  set ρ := testRank A j with hρ
  set c := (univ \ A).orderEmbOfFin hc
  have hcmem : c ⟨ρ, hlt⟩ ∈ univ \ A := (univ \ A).orderEmbOfFin_mem hc _
  have hne : A.orderEmbOfFin h j ≠ c ⟨ρ, hlt⟩ := fun he => by
    have := (mem_sdiff.1 hcmem).2
    rw [← he] at this
    exact this (A.orderEmbOfFin_mem h j)
  rcases lt_or_gt_of_ne hne with h1 | h1
  · exact h1
  · exfalso
    have hsub : insert (c ⟨ρ, hlt⟩) ((univ \ A).filter fun p => p < c ⟨ρ, hlt⟩) ⊆
        (univ \ A).filter fun p => p < A.orderEmbOfFin h j := by
      intro p hp
      rw [mem_insert] at hp
      rw [mem_filter]
      rcases hp with rfl | hp
      · exact ⟨hcmem, h1⟩
      · rw [mem_filter] at hp
        exact ⟨hp.1, hp.2.trans h1⟩
    have hT : ((univ \ A).filter fun p => p < A.orderEmbOfFin h j).card = ρ := by
      rw [hρ, testRank, dif_pos h]
    have := card_le_card hsub
    rw [card_insert_of_notMem (by simp), card_filter_lt_orderEmbOfFin, hT] at this
    simp only at this
    omega

/-- **The order-statistic sandwich, from below**: the calibration unit just under the `j`-th
smallest test unit's calibration rank lies below it. -/
theorem cal_lt_orderEmbOfFin (A : Finset (Fin (n + m))) (h : A.card = m)
    (hc : (univ \ A).card = n) (j : Fin m) (hpos : 0 < testRank A j) (hle : testRank A j ≤ n) :
    (univ \ A).orderEmbOfFin hc ⟨testRank A j - 1, by omega⟩ < A.orderEmbOfFin h j := by
  set ρ := testRank A j with hρ
  set c := (univ \ A).orderEmbOfFin hc
  have hcmem : c ⟨ρ - 1, by omega⟩ ∈ univ \ A := (univ \ A).orderEmbOfFin_mem hc _
  have hne : c ⟨ρ - 1, by omega⟩ ≠ A.orderEmbOfFin h j := fun he => by
    have := (mem_sdiff.1 hcmem).2
    rw [he] at this
    exact this (A.orderEmbOfFin_mem h j)
  rcases lt_or_gt_of_ne hne with h1 | h1
  · exact h1
  · exfalso
    have hsub : ((univ \ A).filter fun p => p < A.orderEmbOfFin h j) ⊆
        (univ \ A).filter fun p => p < c ⟨ρ - 1, by omega⟩ := by
      intro p hp
      rw [mem_filter] at hp ⊢
      exact ⟨hp.1, hp.2.trans h1⟩
    have hT : ((univ \ A).filter fun p => p < A.orderEmbOfFin h j).card = ρ := by
      rw [hρ, testRank, dif_pos h]
    have := card_le_card hsub
    rw [card_filter_lt_orderEmbOfFin, hT] at this
    simp only at this
    omega

/-- The calibration scores in order, with `top` past the last: their `S_(r)`, with
`S_(n+1) = sup s`. -/
noncomputable def calScore (S : Fin (n + m) → ℝ) (top : ℝ) (A : Finset (Fin (n + m)))
    (hc : (univ \ A).card = n) (k : ℕ) : ℝ :=
  if hk : k < n then S ((univ \ A).orderEmbOfFin hc ⟨k, hk⟩) else top

/-- The calibration score just under rank `k`: their `S_(k−1)` in the 1-based count, with
`S_(0) = inf s`, given as `bot`. -/
noncomputable def calScoreBelow (S : Fin (n + m) → ℝ) (bot : ℝ) (A : Finset (Fin (n + m)))
    (hc : (univ \ A).card = n) (k : ℕ) : ℝ :=
  if hk : 0 < k ∧ k ≤ n then
    S ((univ \ A).orderEmbOfFin hc ⟨k - 1, by have := hk.1; have := hk.2; omega⟩)
  else bot

/-- A rank sequence's entries are at most `n`. -/
theorem le_of_mem_rankSeqs {r : Fin m → ℕ} (hr : r ∈ rankSeqs n m) (j : Fin m) : r j ≤ n :=
  Nat.lt_succ_iff.1 (mem_range.1 (Fintype.mem_piFinset.1 (mem_filter.1 hr).1 j))

theorem calScore_mono (S : Fin (n + m) → ℝ) (hS : Monotone S) (top : ℝ) (htop : ∀ p, S p ≤ top)
    (A : Finset (Fin (n + m))) (hc : (univ \ A).card = n) : Monotone (calScore S top A hc) := by
  intro k k' hkk'
  unfold calScore
  by_cases h1 : k < n
  · by_cases h2 : k' < n
    · rw [dif_pos h1, dif_pos h2]
      exact hS (((univ \ A).orderEmbOfFin hc).monotone (Fin.mk_le_mk.2 hkk'))
    · rw [dif_pos h1, dif_neg h2]
      exact htop _
  · have h2 : ¬ k' < n := fun h => h1 (lt_of_le_of_lt hkk' h)
    rw [dif_neg h1, dif_neg h2]

theorem calScoreBelow_mono (S : Fin (n + m) → ℝ) (hS : Monotone S) (bot : ℝ)
    (hbot : ∀ p, bot ≤ S p) (A : Finset (Fin (n + m))) (hc : (univ \ A).card = n) {k k' : ℕ}
    (hkk' : k ≤ k') (hk' : k' ≤ n) :
    calScoreBelow S bot A hc k ≤ calScoreBelow S bot A hc k' := by
  unfold calScoreBelow
  by_cases h1 : 0 < k ∧ k ≤ n
  · have h2 : 0 < k' ∧ k' ≤ n := ⟨lt_of_lt_of_le h1.1 hkk', hk'⟩
    rw [dif_pos h1, dif_pos h2]
    exact hS (((univ \ A).orderEmbOfFin hc).monotone
      (Fin.mk_le_mk.2 (Nat.sub_le_sub_right hkk' 1)))
  · rw [dif_neg h1]
    split_ifs
    · exact hbot _
    · exact le_rfl

/-- **The ranks sandwich the sorted test scores** (the core of their Theorem 1's proof, every
coordinate at once). Units are their positions in score order, `S` monotone between `bot` and
`top`. Ties are broken uniformly at random, which leaves the test set uniform over the
`C(n + m, m)` subsets and moves no calibration order statistic. For any set `I` of rank sequences,
on at least `|H ∩ I|` of the test sets one `r ∈ I` puts the `j`-th smallest test score between
`S_(r_j − 1)` and `S_(r_j)` for every `j`. -/
theorem batch_sandwich (S : Fin (n + m) → ℝ) (hS : Monotone S) (bot top : ℝ)
    (hbot : ∀ p, bot ≤ S p) (htop : ∀ p, S p ≤ top) (I : Finset (Fin m → ℕ)) :
    (rankSeqs n m ∩ I).card ≤ ((powersetCard m (univ : Finset (Fin (n + m)))).filter fun A =>
      ∃ h : A.card = m, ∃ hc : (univ \ A).card = n, ∃ r ∈ rankSeqs n m ∩ I, ∀ j,
        calScoreBelow S bot A hc (r j) ≤ S (A.orderEmbOfFin h j) ∧
        S (A.orderEmbOfFin h j) ≤ calScore S top A hc (r j)).card := by
  rw [← card_testRank_mem I]
  refine card_le_card fun A hA => ?_
  rw [mem_filter, mem_powersetCard] at hA ⊢
  have h := hA.1.2
  have hc : (univ \ A).card = n := by
    rw [card_sdiff_of_subset (subset_univ A), card_univ, Fintype.card_fin, h]
    omega
  refine ⟨hA.1, h, hc, testRank A, mem_inter.2 ⟨testRank_mem_rankSeqs A h, hA.2⟩,
    fun j => ⟨?_, ?_⟩⟩
  · unfold calScoreBelow
    split_ifs with hk
    · exact hS (cal_lt_orderEmbOfFin A h hc j hk.1 hk.2).le
    · exact hbot _
  · unfold calScore
    split_ifs with hlt
    · exact hS (orderEmbOfFin_lt_cal A h hc j hlt).le
    · exact htop _

/-- **Batch predictive inference** (their Theorem 1). The setting is `batch_sandwich`'s, with `h`
monotone in the sorted test scores. For any set `I` of rank sequences, on at least `|H ∩ I|` of the
test sets one `r ∈ I` puts `h(sorted test scores)` between `h(S_(r_1 − 1), …)` and
`h(S_(r_1), …)`. With `I = {q_L ≤ h̃ ≤ q_U}` that is inside their `[B_L, B_U]`, whose min runs over
`{h̃ ≥ q_L}` and max over `{h̃ ≤ q_U}`, both containing `I`. Their `h` on sorted vectors extends to
`v ↦ h(v↑)`, monotone by their Condition 1. -/
theorem batch_interval (S : Fin (n + m) → ℝ) (hS : Monotone S) (bot top : ℝ)
    (hbot : ∀ p, bot ≤ S p) (htop : ∀ p, S p ≤ top) (hh : (Fin m → ℝ) → ℝ) (hmono : Monotone hh)
    (I : Finset (Fin m → ℕ)) :
    (rankSeqs n m ∩ I).card ≤ ((powersetCard m (univ : Finset (Fin (n + m)))).filter fun A =>
      ∃ h : A.card = m, ∃ hc : (univ \ A).card = n, ∃ r ∈ rankSeqs n m ∩ I,
        hh (fun j => calScoreBelow S bot A hc (r j)) ≤ hh (fun j => S (A.orderEmbOfFin h j)) ∧
        hh (fun j => S (A.orderEmbOfFin h j)) ≤ hh (fun j => calScore S top A hc (r j))).card := by
  refine (batch_sandwich S hS bot top hbot htop I).trans (card_le_card fun A hA => ?_)
  rw [mem_filter] at hA ⊢
  obtain ⟨hA, h, hc, r, hr, hs⟩ := hA
  exact ⟨hA, h, hc, r, hr, hmono fun j => (hs j).1, hmono fun j => (hs j).2⟩

/-- From a count of rank sequences to a fraction of test sets: `|H| = C(n + m, m)`. -/
theorem frac_of_card {α : ℝ} {I : Finset (Fin m → ℕ)} {K : ℕ}
    (hI : (1 - α) * (rankSeqs n m).card ≤ (rankSeqs n m ∩ I).card)
    (hK : (rankSeqs n m ∩ I).card ≤ K) : 1 - α ≤ (K : ℝ) / (n + m).choose m := by
  rw [le_div_iff₀ (by exact_mod_cast Nat.choose_pos (Nat.le_add_left m n)),
    ← card_rankSeqs (n := n) (m := m)]
  exact hI.trans (by exact_mod_cast hK)

/-- **Batch coverage** (their Theorem 1): when `I` holds at least `1 − α` of the rank sequences,
as `{q_L ≤ h̃ ≤ q_U}` does by the quantiles' definitions (at most `β` of them below `q_L`, at most
`γ` above `q_U`, `β + γ = α`), the interval holds on at least `1 − α` of the `C(n + m, m)` test
sets. -/
theorem batch_coverage (S : Fin (n + m) → ℝ) (hS : Monotone S) (bot top : ℝ)
    (hbot : ∀ p, bot ≤ S p) (htop : ∀ p, S p ≤ top) (hh : (Fin m → ℝ) → ℝ) (hmono : Monotone hh)
    (I : Finset (Fin m → ℕ)) (α : ℝ)
    (hI : (1 - α) * (rankSeqs n m).card ≤ (rankSeqs n m ∩ I).card) :
    1 - α ≤ (((powersetCard m (univ : Finset (Fin (n + m)))).filter fun A =>
      ∃ h : A.card = m, ∃ hc : (univ \ A).card = n, ∃ r ∈ rankSeqs n m ∩ I,
        hh (fun j => calScoreBelow S bot A hc (r j)) ≤ hh (fun j => S (A.orderEmbOfFin h j)) ∧
        hh (fun j => S (A.orderEmbOfFin h j)) ≤ hh (fun j => calScore S top A hc (r j))).card :
          ℝ) / (n + m).choose m :=
  frac_of_card hI (batch_interval S hS bot top hbot htop hh hmono I)

/-- **Simultaneous sets** (their theorem on joint bounds for several test order statistics). With
`I` the box `w_j ≤ r_j ≤ q_j` (a coordinate left free takes `0` and `n`), on at least `|H ∩ I|` of
the test sets every sorted test score lies in `[S_(w_j − 1), S_(q_j)]`. Their condition
`F_{n,m}(w; q) ≥ (1 − α) C(n + m, m)` is `|H ∩ I| ≥ (1 − α)|H|`, and `frac_of_card` gives
`1 − α`. -/
theorem batch_simultaneous (S : Fin (n + m) → ℝ) (hS : Monotone S) (bot top : ℝ)
    (hbot : ∀ p, bot ≤ S p) (htop : ∀ p, S p ≤ top) (w q : Fin m → ℕ) :
    (rankSeqs n m ∩ (rankSeqs n m).filter fun r => ∀ j, w j ≤ r j ∧ r j ≤ q j).card ≤
      ((powersetCard m (univ : Finset (Fin (n + m)))).filter fun A =>
        ∃ h : A.card = m, ∃ hc : (univ \ A).card = n, ∀ j,
          calScoreBelow S bot A hc (w j) ≤ S (A.orderEmbOfFin h j) ∧
          S (A.orderEmbOfFin h j) ≤ calScore S top A hc (q j)).card := by
  refine (batch_sandwich S hS bot top hbot htop _).trans (card_le_card fun A hA => ?_)
  rw [mem_filter] at hA ⊢
  obtain ⟨hA, h, hc, r, hr, hs⟩ := hA
  obtain ⟨hrH, hrI⟩ := mem_inter.1 hr
  have hb := (mem_filter.1 hrI).2
  refine ⟨hA, h, hc, fun j => ⟨?_, ?_⟩⟩
  · exact (calScoreBelow_mono S hS bot hbot A hc (hb j).1 (le_of_mem_rankSeqs hrH j)).trans
      (hs j).1
  · exact (hs j).2.trans (calScore_mono S hS top htop A hc (hb j).2)

/-- **Exceedances** (their quantile corollary, one-sided, and their selection corollary). On at
least `|H ∩ {r_ζ ≤ k}|` of the test sets, at most `m − 1 − ζ` test scores exceed `S_(k)`, `ζ`
0-based: their `η` false claims at `T̂ = S_(q_η)`, with 1-based `ζ = m − η`. -/
theorem batch_exceed (S : Fin (n + m) → ℝ) (hS : Monotone S) (bot top : ℝ)
    (hbot : ∀ p, bot ≤ S p) (htop : ∀ p, S p ≤ top) (ζ : Fin m) (k : ℕ) :
    (rankSeqs n m ∩ (rankSeqs n m).filter fun r => r ζ ≤ k).card ≤
      ((powersetCard m (univ : Finset (Fin (n + m)))).filter fun A =>
        ∃ h : A.card = m, ∃ hc : (univ \ A).card = n,
          (univ.filter fun j : Fin m => calScore S top A hc k < S (A.orderEmbOfFin h j)).card ≤
            m - 1 - ζ).card := by
  refine (batch_sandwich S hS bot top hbot htop _).trans (card_le_card fun A hA => ?_)
  rw [mem_filter] at hA ⊢
  obtain ⟨hA, h, hc, r, hr, hs⟩ := hA
  have hk := (mem_filter.1 (mem_inter.1 hr).2).2
  refine ⟨hA, h, hc, ?_⟩
  have hsub : (univ.filter fun j : Fin m => calScore S top A hc k < S (A.orderEmbOfFin h j)) ⊆
      Finset.Ioi ζ := by
    intro j hj
    rw [mem_filter] at hj
    rw [Finset.mem_Ioi]
    refine not_le.1 fun hjζ => ?_
    have := (hS ((A.orderEmbOfFin h).monotone hjζ)).trans
      ((hs ζ).2.trans (calScore_mono S hS top htop A hc hk))
    exact absurd hj.2 (not_lt.2 this)
  exact (card_le_card hsub).trans (Fin.card_Ioi ζ).le

/-- Their selection score `s = μ̂ · 1{y ≤ c}`: at a threshold `T ≥ 0` (a score's value, when
`μ̂ ≥ 0`), a false claim, `μ̂ > T` while `y ≤ c`, is exactly a score above `T`. -/
theorem false_claim_iff {μ T y c : ℝ} (hT : 0 ≤ T) :
    (T < μ ∧ y ≤ c) ↔ T < if y ≤ c then μ else 0 := by
  split_ifs with h
  · exact ⟨fun h' => h'.1, fun h' => ⟨h', h⟩⟩
  · exact ⟨fun h' => absurd h'.2 h, fun h' => absurd h' (not_lt.2 hT)⟩

end Batch

/-! ### Slice 18: the guarantee read out

Snell and Griffiths, arXiv 2502.13228, the posterior over the loss; Gressmann et al., arXiv
1801.00753, proper scores. -/

section LossPosterior

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- A sum over `C` is the same after a permutation that fixes everything outside `C`. -/
theorem sum_comp_of_mem_permsOn {C : Finset ι} {σ : Equiv.Perm ι} (hσ : σ ∈ permsOn C)
    (f : ι → ℝ) : ∑ u ∈ C, f (σ u) = ∑ u ∈ C, f u := by
  refine sum_nbij' σ σ.symm (fun u hu => ?_) (fun u hu => ?_) (fun u _ => by simp)
    (fun u _ => by simp) (fun u _ => rfl)
  · exact (mem_iff_of_mem_permsOn hσ u).2 hu
  · exact (mem_iff_of_mem_permsOn hσ (σ.symm u)).1 (by simpa using hu)

/-- The first threshold on a grid meeting a criterion, the last when none does: the `inf` of a
conformal risk control rule. -/
noncomputable def firstFeasible {k : ℕ} (P : Fin (k + 1) → Prop) [DecidablePred P] :
    Fin (k + 1) :=
  if h : (univ.filter P).Nonempty then (univ.filter P).min' h else Fin.last k

theorem firstFeasible_le {k : ℕ} {P : Fin (k + 1) → Prop} [DecidablePred P] {l : Fin (k + 1)}
    (hl : P l) : firstFeasible P ≤ l := by
  have hm : l ∈ univ.filter P := mem_filter.2 ⟨mem_univ _, hl⟩
  unfold firstFeasible
  rw [dif_pos ⟨l, hm⟩]
  exact min'_le _ _ hm

theorem firstFeasible_spec {k : ℕ} {P : Fin (k + 1) → Prop} [DecidablePred P]
    (hlast : P (Fin.last k)) : P (firstFeasible P) := by
  have hne : (univ.filter P).Nonempty := ⟨_, mem_filter.2 ⟨mem_univ _, hlast⟩⟩
  unfold firstFeasible
  rw [dif_pos hne]
  exact (mem_filter.1 (min'_mem _ hne)).2

/-- A weaker criterion is met no later. -/
theorem firstFeasible_mono {k : ℕ} {P R : Fin (k + 1) → Prop} [DecidablePred P] [DecidablePred R]
    (hPR : ∀ l, P l → R l) : firstFeasible R ≤ firstFeasible P := by
  by_cases h : (univ.filter P).Nonempty
  · have hP : P (firstFeasible P) := by
      unfold firstFeasible
      rw [dif_pos h]
      exact (mem_filter.1 (min'_mem _ h)).2
    exact firstFeasible_le (hPR _ hP)
  · have : firstFeasible P = Fin.last k := by
      unfold firstFeasible
      rw [dif_neg h]
    rw [this]
    exact Fin.le_last _

theorem firstFeasible_congr {k : ℕ} {P R : Fin (k + 1) → Prop} [DecidablePred P]
    [DecidablePred R] (hPR : ∀ l, P l ↔ R l) : firstFeasible P = firstFeasible R :=
  le_antisymm (firstFeasible_mono fun l h => (hPR l).2 h)
    (firstFeasible_mono fun l h => (hPR l).1 h)

/-- The conformal risk control rule (Angelopoulos et al.; Snell and Griffiths's `λ_crc`, the
`inf` over `λ` of `E(L⁺) ≤ α`): the first threshold at which the calibration units' losses, with
`B` standing in for the test unit, total at most `α` per unit of `C`. -/
noncomputable def crcRule {Z : Type*} {k : ℕ} (C : Finset ι) (t : ι) (ℓ : Z → Fin (k + 1) → ℝ)
    (B α : ℝ) (z : ι → Z) : Fin (k + 1) :=
  firstFeasible fun l => ∑ i ∈ C.erase t, ℓ (z i) l + B ≤ α * C.card

/-- **Conformal risk control, over orderings** (Snell and Griffiths's Proposition 3.2, through
Angelopoulos et al.'s Theorem 1, arXiv 2208.02814): the posterior mean `E(L⁺)` is calibrated over
replicated calibration sets. Losses nonincreasing in the threshold, at most `B`, and at most `α` at
the last threshold; on a grid, their right-continuity is not needed. Averaged over the orderings of
any data set, the test unit's loss at the rule's threshold is at most `α`. -/
theorem crc_orderings {Z : Type*} {k : ℕ} {C : Finset ι} {t : ι} (ht : t ∈ C)
    (ℓ : Z → Fin (k + 1) → ℝ) (hanti : ∀ v, Antitone (ℓ v)) {B α : ℝ} (hB : ∀ v l, ℓ v l ≤ B)
    (hlast : ∀ v, ℓ v (Fin.last k) ≤ α) (z : ι → Z) :
    ∑ σ ∈ permsOn C, ℓ (z (σ t)) (crcRule C t ℓ B α (z ∘ σ)) ≤ α * (permsOn C).card := by
  set lo := firstFeasible fun l => ∑ u ∈ C, ℓ (z u) l ≤ α * C.card with hlo
  have hlo_spec : ∑ u ∈ C, ℓ (z u) lo ≤ α * C.card := by
    refine firstFeasible_spec (P := fun l => ∑ u ∈ C, ℓ (z u) l ≤ α * C.card) ?_
    calc ∑ u ∈ C, ℓ (z u) (Fin.last k) ≤ ∑ u ∈ C, α := sum_le_sum fun u _ => hlast (z u)
      _ = α * C.card := by rw [sum_const, nsmul_eq_mul, mul_comm]
  have hstep : ∀ σ ∈ permsOn C,
      ℓ (z (σ t)) (crcRule C t ℓ B α (z ∘ σ)) ≤ ℓ (z (σ t)) lo := by
    intro σ hσ
    have hor : (firstFeasible fun l => ∑ u ∈ C, ℓ ((z ∘ σ) u) l ≤ α * C.card) = lo := by
      rw [hlo]
      refine firstFeasible_congr fun l => ?_
      rw [show ∑ u ∈ C, ℓ ((z ∘ σ) u) l = ∑ u ∈ C, ℓ (z u) l from
        sum_comp_of_mem_permsOn hσ (fun u => ℓ (z u) l)]
    have hle : lo ≤ crcRule C t ℓ B α (z ∘ σ) := by
      rw [← hor]
      unfold crcRule
      refine firstFeasible_mono fun l hl => ?_
      beta_reduce at hl ⊢
      rw [← add_sum_erase C _ ht]
      have := hB ((z ∘ σ) t) l
      linarith
    exact hanti _ hle
  have hC : (0 : ℝ) < C.card := by exact_mod_cast card_pos.2 ⟨t, ht⟩
  have hsum : (C.card : ℝ) * ∑ σ ∈ permsOn C, ℓ (z (σ t)) lo =
      (permsOn C).card * ∑ u ∈ C, ℓ (z u) lo := card_mul_sum_permsOn ht (fun u => ℓ (z u) lo)
  have h1 := sum_le_sum hstep
  have h2 : (C.card : ℝ) * ∑ σ ∈ permsOn C, ℓ (z (σ t)) lo ≤
      (C.card : ℝ) * (α * (permsOn C).card) := by
    rw [hsum]
    calc ((permsOn C).card : ℝ) * ∑ u ∈ C, ℓ (z u) lo ≤ (permsOn C).card * (α * C.card) :=
          mul_le_mul_of_nonneg_left hlo_spec (Nat.cast_nonneg _)
      _ = C.card * (α * (permsOn C).card) := by ring
  have h3 := le_of_mul_le_mul_left h2 hC
  linarith

omit [Fintype ι] [DecidableEq ι] in
/-- **The `β`-bound's calibration at one calibration unit** (Gaffke's bound, which is Snell and
Griffiths's `b*_β`). With one loss `ℓ` their posterior is `L⁺ = ℓ + U (B − ℓ)`, `U` uniform, so
`b*_β = ℓ + β (B − ℓ)`. For a loss law on finitely many values at most `B`, the bound is at least
the mean loss `J` with chance at least `β`. Every number of units is `gaffke_calibration`
(Vlassis and Thomas's Theorem 1, arXiv 2607.08415), through `P(L⁺ ≤ J) = 1 − K(X)` with
`X_i = (B − ℓ_i)/(B − J)`. -/
theorem gaffke_one {κ : Type*} [Fintype κ] (p v : κ → ℝ) (hp : ∀ k, 0 ≤ p k)
    (hp1 : ∑ k, p k = 1) {B β : ℝ} (hv : ∀ k, v k ≤ B) (hβ : β ≤ 1) :
    β ≤ ∑ k, p k * if (∑ k', p k' * v k') ≤ v k + β * (B - v k) then 1 else 0 := by
  classical
  set J := ∑ k', p k' * v k' with hJ
  set G := univ.filter fun k => J ≤ v k + β * (B - v k) with hG
  set F := univ.filter fun k => ¬ (J ≤ v k + β * (B - v k)) with hF
  have hsplit : ∀ f : κ → ℝ, ∑ k, f k = ∑ k ∈ G, f k + ∑ k ∈ F, f k :=
    fun f => (sum_filter_add_sum_filter_not univ _ f).symm
  have hrhs : ∑ k, p k * (if J ≤ v k + β * (B - v k) then 1 else 0) = 1 - ∑ k ∈ F, p k := by
    have h1 := hsplit p
    rw [hp1] at h1
    simp only [mul_ite, mul_one, mul_zero]
    rw [← sum_filter]
    linarith
  rw [hrhs]
  set P := ∑ k ∈ F, p k with hP
  have hJB : J ≤ B :=
    calc J = ∑ k, p k * v k := rfl
      _ ≤ ∑ k, p k * B := sum_le_sum fun k _ => mul_le_mul_of_nonneg_left (hv k) (hp k)
      _ = B := by rw [← sum_mul, hp1, one_mul]
  have hb : J ≤ ∑ k ∈ F, p k * v k + B * (1 - P) := by
    have h1 := hsplit fun k => p k * v k
    have h2 := hsplit p
    rw [hp1] at h2
    have h3 : ∑ k ∈ G, p k * v k ≤ ∑ k ∈ G, p k * B :=
      sum_le_sum fun k _ => mul_le_mul_of_nonneg_left (hv k) (hp k)
    rw [← sum_mul, show ∑ k ∈ G, p k = 1 - P by linarith] at h3
    have : J = ∑ k ∈ G, p k * v k + ∑ k ∈ F, p k * v k := h1
    nlinarith
  by_cases hP0 : P = 0
  · rw [hP0]
    linarith
  obtain ⟨k, hkF, hpk⟩ : ∃ k ∈ F, 0 < p k := by
    by_contra hc
    exact hP0 (sum_eq_zero fun k hk => le_antisymm (not_lt.1 fun h => hc ⟨k, hk, h⟩) (hp k))
  have hstrict : ∑ k ∈ F, p k * (v k + β * (B - v k)) < ∑ k ∈ F, p k * J :=
    sum_lt_sum (fun k hk => mul_le_mul_of_nonneg_left (not_le.1 (mem_filter.1 hk).2).le (hp k))
      ⟨k, hkF, mul_lt_mul_of_pos_left (not_le.1 (mem_filter.1 hkF).2) hpk⟩
  have hexp : ∑ k ∈ F, p k * (v k + β * (B - v k)) =
      (1 - β) * ∑ k ∈ F, p k * v k + β * B * P := by
    rw [hP, mul_sum, mul_sum, ← sum_add_distrib]
    exact sum_congr rfl fun k _ => by ring
  rw [hexp, ← sum_mul] at hstrict
  have hmul := mul_le_mul_of_nonneg_left hb (sub_nonneg.2 hβ)
  by_contra hc
  have hc' : 1 - P < β := not_le.1 hc
  nlinarith [mul_nonneg (sub_nonneg.2 hJB) (sub_nonneg.2 hc'.le)]

/-- **A data-chosen threshold fails only where a fixed one does** (the step behind Snell and
Griffiths's `λ^β_hpd` controlling the risk). With the bound `b` nonincreasing in the threshold,
`l` chosen with `b l ≤ α`, and `l₀` the last threshold whose risk `J` exceeds `α`: risk above `α`
at `l` puts the bound below the risk at `l₀`, an event of chance at most `1 − β` by the fixed
threshold's calibration (`gaffke_miss_le`); `highProbability_risk` puts the two together. -/
theorem fail_le_fixed {Λ : Type*} [Preorder Λ] {J b : Λ → ℝ} (hb : Antitone b) {α : ℝ}
    {l l₀ : Λ} (hl : b l ≤ α) (hJ₀ : α < J l₀) (hl₀ : ∀ l', α < J l' → l' ≤ l₀)
    (hfail : α < J l) : b l₀ < J l₀ :=
  lt_of_le_of_lt ((hb (hl₀ l hfail)).trans hl) hJ₀

end LossPosterior

/-! ### Slice 18's follow-up: the calibration at `n` units

Vlassis and Thomas, arXiv 2607.08415, prove Gaffke's conjecture: independent nonnegative `X_i`
with means at most `1` have `P(K(X) ≤ α) ≤ α`, `K(x) = P(∑ x_i D_i ≤ 1)` with `D ~ Dir(1, …, 1)`
over the `x_i` and an added zero. With `X_i = (B − ℓ_i)/(B − J)` it is the `β`-bound's calibration
at every number of units, and with `fail_le_fixed` it is what `RiskControl.HighProbability`
reports.

Laws on the line built by convolution carry the exponentials: `sumLaw m` is the law of
`∑_{r ∈ m} r E_r`, and `stat m = P(∑ r E_r ≥ 0)`. The race identity (`stat_race`) is the one
analytic computation. The same-sign and flux identities follow from it by algebra. The
likelihood-ratio inequality of their exponential transfer lemma (`lrAbove_edge`) follows by adding
one exponential at a time, memorylessness splitting each integral (`smooth_ordered`), so no
log-concavity theorem is needed. The chain-measure induction is `local_insertion`,
`chain_dominates` and `twoPoint_valid`. The reduction to finite laws is `twoPoint_mixture`,
`sum_prod_twoPoint` and `gaffkeK_valid`. The calibration is `gaffke_miss_le`, `gaffke_calibration`
and `highProbability_risk`. -/

section VlassisThomas

open MeasureTheory Set Real Filter
open scoped ENNReal

/-- The unit exponential's density. -/
noncomputable def unitExpDensity (s : ℝ) : ℝ≥0∞ := ENNReal.ofReal (if 0 ≤ s then exp (-s) else 0)

theorem measurable_unitExpDensity : Measurable unitExpDensity :=
  ENNReal.measurable_ofReal.comp
    (Measurable.ite (measurableSet_Ici : MeasurableSet (Set.Ici (0 : ℝ))) (by fun_prop)
      measurable_const)

/-- The unit exponential law, written by its density. -/
noncomputable def unitExp : Measure ℝ := volume.withDensity unitExpDensity

theorem lintegral_unitExp {f : ℝ → ℝ≥0∞} (hf : Measurable f) :
    ∫⁻ s, f s ∂unitExp = ∫⁻ s, unitExpDensity s * f s := by
  rw [unitExp, lintegral_withDensity_eq_lintegral_mul _ measurable_unitExpDensity hf]
  rfl

theorem ae_unitExp_nonneg : ∀ᵐ s ∂unitExp, 0 ≤ s := by
  rw [unitExp, ae_withDensity_iff measurable_unitExpDensity]
  refine Eventually.of_forall fun s hs => ?_
  by_contra h
  simp [unitExpDensity, h] at hs

/-- The unit exponential's Laplace transform. -/
theorem lintegral_unitExp_exp {l : ℝ} (hl : 0 ≤ l) :
    ∫⁻ s, ENNReal.ofReal (exp (-(l * s))) ∂unitExp = ENNReal.ofReal (1 / (1 + l)) := by
  have hm : Measurable fun s : ℝ => ENNReal.ofReal (exp (-(l * s))) := by fun_prop
  rw [lintegral_unitExp hm]
  have h1 : ∀ s, unitExpDensity s * ENNReal.ofReal (exp (-(l * s))) =
      (Set.Ioi 0).indicator (fun s => ENNReal.ofReal (exp ((-(1 + l)) * s))) s +
        ({0} : Set ℝ).indicator (fun _ => 1) s := by
    intro s
    rcases lt_trichotomy s 0 with hs | rfl | hs
    · simp [unitExpDensity, hs.not_ge, hs.ne,
        Set.indicator_of_notMem (show s ∉ Set.Ioi (0 : ℝ) from fun h => by simp at h; linarith)]
    · simp [unitExpDensity]
    · rw [Set.indicator_of_mem (show s ∈ Set.Ioi (0 : ℝ) from hs),
        Set.indicator_of_notMem (show s ∉ ({0} : Set ℝ) from fun h => by simp at h; linarith)]
      simp only [unitExpDensity, if_pos hs.le, add_zero]
      rw [← ENNReal.ofReal_mul (exp_pos _).le, ← exp_add]
      congr 2
      ring
  simp_rw [h1]
  have hmeas :
      Measurable ((Set.Ioi (0 : ℝ)).indicator fun s => ENNReal.ofReal (exp ((-(1 + l)) * s))) :=
    (by fun_prop : Measurable fun s : ℝ => ENNReal.ofReal (exp ((-(1 + l)) * s))).indicator
      measurableSet_Ioi
  rw [lintegral_add_left hmeas,
    lintegral_indicator measurableSet_Ioi, lintegral_indicator_const (measurableSet_singleton 0),
    Real.volume_singleton, mul_zero, add_zero,
    ← ofReal_integral_eq_lintegral_ofReal (integrableOn_exp_mul_Ioi (by linarith) 0)
      (Eventually.of_forall fun s => (exp_pos _).le),
    integral_exp_mul_Ioi (by linarith) 0]
  congr 1
  simp only [mul_zero, exp_zero]
  field_simp

instance : IsProbabilityMeasure unitExp := by
  refine ⟨?_⟩
  have h := lintegral_unitExp_exp (l := 0) le_rfl
  simp only [zero_mul, neg_zero, exp_zero, ENNReal.ofReal_one, add_zero, div_one] at h
  rw [← lintegral_one]
  exact h

/-- **Memorylessness.** Past `δ ≥ 0`, the unit exponential is `δ` plus a fresh one, with chance
`e^{−δ}`. -/
theorem lintegral_unitExp_Ici {f : ℝ → ℝ≥0∞} (hf : Measurable f) {δ : ℝ} (hδ : 0 ≤ δ) :
    ∫⁻ s, (Set.Ici δ).indicator f s ∂unitExp =
      ENNReal.ofReal (exp (-δ)) * ∫⁻ s, f (s + δ) ∂unitExp := by
  have hfδ : Measurable fun s => f (s + δ) := hf.comp (measurable_add_const δ)
  rw [lintegral_unitExp (hf.indicator measurableSet_Ici), lintegral_unitExp hfδ,
    ← lintegral_const_mul' _ _ ENNReal.ofReal_ne_top,
    ← lintegral_add_right_eq_self (fun s => unitExpDensity s * (Set.Ici δ).indicator f s) δ]
  refine lintegral_congr fun s => ?_
  by_cases hs : 0 ≤ s
  · have h1 : s + δ ∈ Set.Ici δ := by simp [hs]
    rw [Set.indicator_of_mem h1]
    simp only [unitExpDensity, if_pos hs, if_pos (add_nonneg hs hδ)]
    rw [← mul_assoc, ← ENNReal.ofReal_mul (exp_pos _).le, ← exp_add]
    congr 3
    ring
  · have h1 : s + δ ∉ Set.Ici δ := by simp only [Set.mem_Ici]; linarith
    rw [Set.indicator_of_notMem h1]
    simp [unitExpDensity, hs]

theorem unitExp_Ici {δ : ℝ} (hδ : 0 ≤ δ) : unitExp (Set.Ici δ) = ENNReal.ofReal (exp (-δ)) := by
  have h := lintegral_unitExp_Ici (f := fun _ => (1 : ℝ≥0∞)) measurable_const hδ
  rw [lintegral_indicator_const measurableSet_Ici, one_mul] at h
  simpa using h

theorem unitExp_Ici_of_nonpos {δ : ℝ} (hδ : δ ≤ 0) : unitExp (Set.Ici δ) = 1 := by
  refine le_antisymm prob_le_one ?_
  calc (1 : ℝ≥0∞) = unitExp (Set.Ici 0) := by rw [unitExp_Ici le_rfl]; simp
    _ ≤ unitExp (Set.Ici δ) := measure_mono (Set.Ici_subset_Ici.2 hδ)

theorem unitExp_Iio_add {δ : ℝ} (hδ : 0 ≤ δ) :
    unitExp (Set.Iio δ) + ENNReal.ofReal (exp (-δ)) = 1 := by
  rw [← unitExp_Ici hδ, ← Set.compl_Ici, measure_compl measurableSet_Ici (measure_ne_top _ _),
    measure_univ]
  exact tsub_add_cancel_of_le prob_le_one

theorem unitExp_Iio_of_nonpos {δ : ℝ} (hδ : δ ≤ 0) : unitExp (Set.Iio δ) = 0 := by
  rw [← Set.compl_Ici, prob_compl_eq_one_sub measurableSet_Ici, unitExp_Ici_of_nonpos hδ, tsub_self]

theorem unitExp_singleton (x : ℝ) : unitExp {x} = 0 := by
  rw [unitExp]
  exact measure_singleton x

/-- The law of `r E`, `E` a unit exponential. -/
noncomputable def expLaw (r : ℝ) : Measure ℝ := unitExp.map (fun s => r * s)

instance (r : ℝ) : IsProbabilityMeasure (expLaw r) :=
  Measure.isProbabilityMeasure_map (measurable_const_mul r).aemeasurable

theorem lintegral_expLaw (r : ℝ) {f : ℝ → ℝ≥0∞} (hf : Measurable f) :
    ∫⁻ x, f x ∂(expLaw r) = ∫⁻ s, f (r * s) ∂unitExp :=
  lintegral_map hf (measurable_const_mul r)

theorem expLaw_apply (r : ℝ) {S : Set ℝ} (hS : MeasurableSet S) :
    expLaw r S = unitExp ((fun s => r * s) ⁻¹' S) :=
  Measure.map_apply (measurable_const_mul r) hS

theorem conv_Ici (μ ν : Measure ℝ) [SFinite ν] (y : ℝ) :
    (μ ∗ ν) (Set.Ici y) = ∫⁻ x, ν (Set.Ici (y - x)) ∂μ := by
  rw [← lintegral_indicator_one measurableSet_Ici,
    Measure.lintegral_conv (measurable_one.indicator measurableSet_Ici)]
  refine lintegral_congr fun x => ?_
  rw [← lintegral_indicator_one measurableSet_Ici]
  refine lintegral_congr fun z => ?_
  simp only [Set.indicator_apply, Set.mem_Ici, Pi.one_apply]
  by_cases h : y ≤ x + z
  · rw [if_pos h, if_pos (by linarith)]
  · rw [if_neg h, if_neg (by intro h'; exact h (by linarith))]

/-- The tail of `a E`, `a > 0`. -/
theorem expLaw_Ici_of_pos {a : ℝ} (ha : 0 < a) {z : ℝ} (hz : 0 ≤ z) :
    expLaw a (Set.Ici z) = ENNReal.ofReal (exp (-(z / a))) := by
  rw [expLaw_apply a measurableSet_Ici,
    show (fun s => a * s) ⁻¹' Set.Ici z = Set.Ici (z / a) by
      ext s; simp [div_le_iff₀ ha, mul_comm]]
  exact unitExp_Ici (div_nonneg hz ha.le)

theorem expLaw_Ici_of_pos_of_nonpos {a : ℝ} (ha : 0 < a) {z : ℝ} (hz : z ≤ 0) :
    expLaw a (Set.Ici z) = 1 := by
  rw [expLaw_apply a measurableSet_Ici,
    show (fun s => a * s) ⁻¹' Set.Ici z = Set.Ici (z / a) by
      ext s; simp [div_le_iff₀ ha, mul_comm]]
  exact unitExp_Ici_of_nonpos (div_nonpos_of_nonpos_of_nonneg hz ha.le)

/-- The tail of `−b E`, `b > 0`. -/
theorem expLaw_neg_Ici_add {b : ℝ} (hb : 0 < b) {z : ℝ} (hz : z ≤ 0) :
    expLaw (-b) (Set.Ici z) + ENNReal.ofReal (exp (z / b)) = 1 := by
  have hpre : (fun s => -b * s) ⁻¹' Set.Ici z = Set.Iic (-z / b) := by
    ext s
    simp only [Set.mem_preimage, Set.mem_Ici, Set.mem_Iic, le_div_iff₀ hb]
    constructor <;> intro h <;> linarith
  rw [expLaw_apply _ measurableSet_Ici, hpre, ← Set.Iio_union_right,
    measure_union (by simp) (measurableSet_singleton _), unitExp_singleton, add_zero,
    show z / b = -(-z / b) by ring]
  exact unitExp_Iio_add (div_nonneg (neg_nonneg.2 hz) hb.le)

theorem expLaw_neg_Ici_of_nonneg {b : ℝ} (hb : 0 < b) {z : ℝ} (hz : 0 ≤ z) :
    expLaw (-b) (Set.Ici z) = 0 := by
  have hpre : (fun s => -b * s) ⁻¹' Set.Ici z ⊆ Set.Iio 0 ∪ {0} := by
    intro s hs
    simp only [Set.mem_preimage, Set.mem_Ici] at hs
    rcases lt_or_ge s 0 with h | h
    · exact Or.inl h
    · exact Or.inr (by simp only [Set.mem_singleton_iff]; nlinarith)
  rw [expLaw_apply _ measurableSet_Ici]
  exact measure_mono_null hpre (measure_union_null (unitExp_Iio_of_nonpos le_rfl)
    (unitExp_singleton 0))

theorem expLaw_zero : expLaw 0 = Measure.dirac 0 := by
  rw [expLaw, show (fun s : ℝ => 0 * s) = fun _ => 0 from funext fun s => zero_mul s,
    Measure.map_const, measure_univ, one_smul]

theorem expLaw_Ici_mono {r r' : ℝ} (h : r ≤ r') (z : ℝ) :
    expLaw r (Set.Ici z) ≤ expLaw r' (Set.Ici z) := by
  rw [expLaw_apply r measurableSet_Ici, expLaw_apply r' measurableSet_Ici]
  refine measure_mono_ae (ae_unitExp_nonneg.mono fun s hs hz => ?_)
  change z ≤ r * s at hz
  change z ≤ r' * s
  have := mul_le_mul_of_nonneg_right h hs
  linarith

theorem expLaw_singleton {r : ℝ} (hr : r ≠ 0) (z : ℝ) : expLaw r {z} = 0 := by
  rw [expLaw_apply r (measurableSet_singleton z)]
  refine measure_mono_null (fun s hs => ?_) (unitExp_singleton (z / r))
  simp only [Set.mem_preimage, Set.mem_singleton_iff] at hs ⊢
  rw [← hs]
  field_simp

/-- **The race of `a E` and `b E'`.** `aE − bE'` is `aE''` with chance `a/(a + b)` and `−bE''`
otherwise; on the half-lines `[y, ∞)`, scaled by `a + b`. -/
theorem race_Ici {a b : ℝ} (ha : 0 < a) (hb : 0 < b) (y : ℝ) :
    ENNReal.ofReal (a + b) * (expLaw (-b) ∗ expLaw a) (Set.Ici y) =
      ENNReal.ofReal a * expLaw a (Set.Ici y) + ENNReal.ofReal b * expLaw (-b) (Set.Ici y) := by
  have hmono : Monotone fun u => expLaw a (Set.Ici (y - u)) :=
    fun u v huv => measure_mono (Set.Ici_subset_Ici.2 (by linarith))
  rw [conv_Ici, lintegral_expLaw _ hmono.measurable]
  have hlap := lintegral_unitExp_exp (div_nonneg hb.le ha.le)
  have hfrac : ENNReal.ofReal (a + b) * ENNReal.ofReal (1 / (1 + b / a)) = ENNReal.ofReal a := by
    rw [← ENNReal.ofReal_mul (by positivity)]
    congr 1
    field_simp
  rcases le_or_gt 0 y with hy | hy
  · have hpt : ∀ᵐ t ∂unitExp, expLaw a (Set.Ici (y - -b * t)) =
        ENNReal.ofReal (exp (-(y / a))) * ENNReal.ofReal (exp (-(b / a * t))) :=
      ae_unitExp_nonneg.mono fun t ht => by
        rw [expLaw_Ici_of_pos ha (by nlinarith), ← ENNReal.ofReal_mul (exp_pos _).le, ← exp_add]
        congr 2
        field_simp
        ring
    rw [lintegral_congr_ae hpt, lintegral_const_mul' _ _ ENNReal.ofReal_ne_top, hlap,
      expLaw_neg_Ici_of_nonneg hb hy, mul_zero, add_zero, expLaw_Ici_of_pos ha hy, mul_comm,
      mul_assoc, mul_comm (ENNReal.ofReal (1 / (1 + b / a))), hfrac, mul_comm]
  · set δ := -y / b with hδ
    have hδ0 : 0 < δ := div_pos (neg_pos.2 hy) hb
    have hsplit : (fun t => expLaw a (Set.Ici (y - -b * t))) = fun t =>
        (Set.Iio δ).indicator (fun _ => (1 : ℝ≥0∞)) t +
          (Set.Ici δ).indicator (fun t => ENNReal.ofReal (exp (-((y + b * t) / a)))) t := by
      funext t
      rcases lt_or_ge t δ with ht | ht
      · rw [Set.indicator_of_mem (show t ∈ Set.Iio δ from ht),
          Set.indicator_of_notMem (show t ∉ Set.Ici δ from fun h => (not_le.2 ht) h), add_zero]
        refine expLaw_Ici_of_pos_of_nonpos ha ?_
        have : b * t < -y := by rw [hδ, lt_div_iff₀ hb] at ht; linarith
        linarith
      · rw [Set.indicator_of_notMem (show t ∉ Set.Iio δ from not_lt.2 ht), Set.indicator_of_mem
          (show t ∈ Set.Ici δ from ht), zero_add, expLaw_Ici_of_pos ha]
        · congr 3; ring
        · have : -y ≤ b * t := by rw [hδ, div_le_iff₀ hb] at ht; linarith
          linarith
    have hm2 : Measurable fun t : ℝ => ENNReal.ofReal (exp (-((y + b * t) / a))) := by fun_prop
    rw [hsplit, lintegral_add_left (measurable_const.indicator measurableSet_Iio),
      lintegral_indicator_const measurableSet_Iio, one_mul,
      lintegral_unitExp_Ici hm2 hδ0.le]
    have hshift : (fun t => ENNReal.ofReal (exp (-((y + b * (t + δ)) / a)))) =
        fun t => ENNReal.ofReal (exp (-(b / a * t))) := by
      funext t
      congr 3
      rw [hδ]
      field_simp
      ring
    rw [hshift, hlap]
    have hE : ENNReal.ofReal (exp (-δ)) = ENNReal.ofReal (exp (y / b)) := by
      rw [hδ]; congr 2; ring
    have hP := unitExp_Iio_add hδ0.le
    have hQ := expLaw_neg_Ici_add hb hy.le
    rw [← hE] at hQ
    have hPQ : unitExp (Set.Iio δ) = expLaw (-b) (Set.Ici y) :=
      (ENNReal.add_left_inj ENNReal.ofReal_ne_top).1 (hP.trans hQ.symm)
    rw [expLaw_Ici_of_pos_of_nonpos ha hy.le, ← hPQ]
    calc ENNReal.ofReal (a + b) *
          (unitExp (Set.Iio δ) + ENNReal.ofReal (exp (-δ)) * ENNReal.ofReal (1 / (1 + b / a)))
        = (ENNReal.ofReal a + ENNReal.ofReal b) * unitExp (Set.Iio δ) + ENNReal.ofReal (exp (-δ)) *
            (ENNReal.ofReal (a + b) * ENNReal.ofReal (1 / (1 + b / a))) := by
          rw [ENNReal.ofReal_add ha.le hb.le]; ring
      _ = ENNReal.ofReal a * (unitExp (Set.Iio δ) + ENNReal.ofReal (exp (-δ))) +
            ENNReal.ofReal b * unitExp (Set.Iio δ) := by rw [hfrac]; ring
      _ = _ := by rw [hP, mul_one]

/-- The adding of an independent `r E`. -/
noncomputable def addExp (r : ℝ) (μ : Measure ℝ) : Measure ℝ := μ ∗ expLaw r

instance : LeftCommutative addExp :=
  ⟨fun a b μ => by
    simp only [addExp]
    rw [Measure.conv_assoc, Measure.conv_assoc, Measure.conv_comm (expLaw b)]⟩

/-- The law of `∑_{r ∈ m} r E_r`, the `E_r` independent unit exponentials. -/
noncomputable def sumLaw (m : Multiset ℝ) : Measure ℝ := Multiset.foldr addExp (Measure.dirac 0) m

theorem sumLaw_zero : sumLaw 0 = Measure.dirac 0 := Multiset.foldr_zero _ _

theorem sumLaw_cons (r : ℝ) (m : Multiset ℝ) : sumLaw (r ::ₘ m) = sumLaw m ∗ expLaw r :=
  Multiset.foldr_cons _ _ _ _

instance (m : Multiset ℝ) : IsProbabilityMeasure (sumLaw m) := by
  induction m using Multiset.induction_on with
  | empty => rw [sumLaw_zero]; infer_instance
  | cons r m ih => rw [sumLaw_cons]; infer_instance

/-- `P(∑_{r ∈ m} r E_r ≥ 0)`. -/
noncomputable def stat (m : Multiset ℝ) : ℝ := (sumLaw m).real (Set.Ici 0)

theorem stat_le_one (m : Multiset ℝ) : stat m ≤ 1 := measureReal_le_one

theorem stat_nonneg (m : Multiset ℝ) : 0 ≤ stat m := measureReal_nonneg

theorem sumLaw_cons_Ici (r : ℝ) (m : Multiset ℝ) :
    sumLaw (r ::ₘ m) (Set.Ici 0) = ∫⁻ x, expLaw r (Set.Ici (-x)) ∂(sumLaw m) := by
  rw [sumLaw_cons, conv_Ici]
  simp

theorem sumLaw_cons_cons_Ici (r r' : ℝ) (m : Multiset ℝ) :
    sumLaw (r ::ₘ r' ::ₘ m) (Set.Ici 0) =
      ∫⁻ x, (expLaw r' ∗ expLaw r) (Set.Ici (-x)) ∂(sumLaw m) := by
  rw [sumLaw_cons, sumLaw_cons, Measure.conv_assoc, conv_Ici]
  simp

theorem measurable_expLaw_Ici_neg (r : ℝ) : Measurable fun x => expLaw r (Set.Ici (-x)) :=
  (show Monotone fun x => expLaw r (Set.Ici (-x)) from
    fun _ _ huv => measure_mono (Set.Ici_subset_Ici.2 (neg_le_neg huv))).measurable

/-- **The race identity** (`(a + b) P(aE − bE' + Y ≥ 0) = a P(aE + Y ≥ 0) + b P(Y − bE ≥ 0)`),
the one analytic computation the induction needs. -/
theorem stat_race {a b : ℝ} (ha : 0 < a) (hb : 0 < b) (Y : Multiset ℝ) :
    (a + b) * stat (a ::ₘ (-b) ::ₘ Y) = a * stat (a ::ₘ Y) + b * stat ((-b) ::ₘ Y) := by
  have h : ENNReal.ofReal (a + b) * sumLaw (a ::ₘ (-b) ::ₘ Y) (Set.Ici 0) =
      ENNReal.ofReal a * sumLaw (a ::ₘ Y) (Set.Ici 0) +
        ENNReal.ofReal b * sumLaw ((-b) ::ₘ Y) (Set.Ici 0) := by
    rw [sumLaw_cons_cons_Ici, sumLaw_cons_Ici, sumLaw_cons_Ici,
      ← lintegral_const_mul' _ _ ENNReal.ofReal_ne_top,
      ← lintegral_const_mul' _ _ ENNReal.ofReal_ne_top,
      ← lintegral_const_mul' _ _ ENNReal.ofReal_ne_top,
      ← lintegral_add_left ((measurable_expLaw_Ici_neg a).const_mul _)]
    exact lintegral_congr fun x => race_Ici ha hb (-x)
  have h2 := congrArg ENNReal.toReal h
  rw [ENNReal.toReal_mul, ENNReal.toReal_add (ENNReal.mul_ne_top ENNReal.ofReal_ne_top
      (measure_ne_top _ _)) (ENNReal.mul_ne_top ENNReal.ofReal_ne_top (measure_ne_top _ _)),
    ENNReal.toReal_mul, ENNReal.toReal_mul, ENNReal.toReal_ofReal (by positivity),
    ENNReal.toReal_ofReal ha.le, ENNReal.toReal_ofReal hb.le] at h2
  exact h2

/-- Two positive coefficients (from the race identity, with an auxiliary `−E`). -/
theorem stat_same_pos {a b : ℝ} (ha : 0 < a) (hb : 0 < b) (Y : Multiset ℝ) :
    (a - b) * stat (a ::ₘ b ::ₘ Y) = a * stat (a ::ₘ Y) - b * stat (b ::ₘ Y) := by
  have e1 := stat_race ha one_pos (b ::ₘ Y)
  have e2 := stat_race hb one_pos (a ::ₘ Y)
  have e3 := stat_race ha one_pos Y
  have e4 := stat_race hb one_pos Y
  rw [Multiset.cons_swap (-1 : ℝ) b] at e1
  rw [Multiset.cons_swap (-1 : ℝ) a, Multiset.cons_swap b a, Multiset.cons_swap b a Y] at e2
  linear_combination (-(b + 1)) * e1 + (a + 1) * e2 + e3 - e4

/-- Two negative coefficients (from the race identity, with an auxiliary `E`). -/
theorem stat_same_neg {a b : ℝ} (ha : 0 < a) (hb : 0 < b) (Y : Multiset ℝ) :
    (a - b) * stat ((-a) ::ₘ (-b) ::ₘ Y) = a * stat ((-a) ::ₘ Y) - b * stat ((-b) ::ₘ Y) := by
  have e1 := stat_race one_pos ha ((-b) ::ₘ Y)
  have e2 := stat_race one_pos hb ((-a) ::ₘ Y)
  have e3 := stat_race one_pos ha Y
  have e4 := stat_race one_pos hb Y
  rw [Multiset.cons_swap (-b) (-a)] at e2
  linear_combination (-(1 + b)) * e1 + (1 + a) * e2 + e3 - e4

theorem stat_cons_le {r r' : ℝ} (h : r ≤ r') (m : Multiset ℝ) :
    stat (r ::ₘ m) ≤ stat (r' ::ₘ m) := by
  unfold stat
  rw [measureReal_def, measureReal_def]
  refine ENNReal.toReal_mono (measure_ne_top _ _) ?_
  rw [sumLaw_cons_Ici, sumLaw_cons_Ici]
  exact lintegral_mono fun x => expLaw_Ici_mono h (-x)

/-- `stat` rises with every coefficient. -/
theorem stat_map_mono {α : Type*} (base : Multiset ℝ) (m : Multiset α) {f g : α → ℝ}
    (hfg : ∀ x ∈ m, f x ≤ g x) : stat (base + m.map f) ≤ stat (base + m.map g) := by
  induction m using Multiset.induction_on generalizing base with
  | empty => simp
  | cons x m ih =>
    simp only [Multiset.map_cons, Multiset.add_cons]
    calc stat (f x ::ₘ (base + m.map f)) ≤ stat (g x ::ₘ (base + m.map f)) :=
          stat_cons_le (hfg x (Multiset.mem_cons_self _ _)) _
      _ = stat ((g x ::ₘ base) + m.map f) := by rw [Multiset.cons_add]
      _ ≤ stat ((g x ::ₘ base) + m.map g) := ih _ fun y hy => hfg y (Multiset.mem_cons_of_mem hy)
      _ = stat (g x ::ₘ (base + m.map g)) := by rw [Multiset.cons_add]

theorem stat_cons_zero (m : Multiset ℝ) : stat (0 ::ₘ m) = stat m := by
  unfold stat
  rw [sumLaw_cons, expLaw_zero, Measure.conv_dirac_zero]

theorem stat_add_zeros {α : Type*} (base : Multiset ℝ) (m : Multiset α) :
    stat (base + m.map fun _ => (0 : ℝ)) = stat base := by
  induction m using Multiset.induction_on with
  | empty => simp
  | cons x m ih => rw [Multiset.map_cons, Multiset.add_cons, stat_cons_zero, ih]

theorem stat_zero : stat 0 = 1 := by
  rw [stat, sumLaw_zero, measureReal_def, Measure.dirac_apply_of_mem (Set.mem_Ici.2 le_rfl),
    ENNReal.toReal_one]

/-- Nonnegative coefficients: the sum is nonnegative surely. -/
theorem stat_eq_one_of_nonneg {m : Multiset ℝ} (hm : ∀ r ∈ m, 0 ≤ r) : stat m = 1 := by
  refine le_antisymm (stat_le_one m) ?_
  have h := stat_map_mono 0 m (f := fun _ => (0 : ℝ)) (g := id) hm
  rwa [stat_add_zeros, stat_zero, zero_add, Multiset.map_id] at h

/-- A negative coefficient and no positive one: the sum is negative surely. -/
theorem stat_eq_zero_of_nonpos {r : ℝ} (hr : r < 0) {m : Multiset ℝ} (hm : ∀ x ∈ m, x ≤ 0) :
    stat (r ::ₘ m) = 0 := by
  refine le_antisymm ?_ (stat_nonneg _)
  have h := stat_map_mono {r} m (f := id) (g := fun _ => (0 : ℝ)) hm
  rw [Multiset.map_id, stat_add_zeros, Multiset.singleton_add] at h
  refine h.trans (le_of_eq ?_)
  have hb : 0 < -r := neg_pos.2 hr
  unfold stat
  rw [← Multiset.cons_zero, sumLaw_cons, sumLaw_zero, Measure.dirac_zero_conv, measureReal_def,
    show r = -(-r) by ring, expLaw_neg_Ici_of_nonneg hb le_rfl, ENNReal.toReal_zero]

/-- The window above zero: `e^{−x/d}` on `x ≥ 0`. -/
noncomputable def winHi (d x : ℝ) : ℝ≥0∞ :=
  (Set.Ici 0).indicator (fun x => ENNReal.ofReal (exp (-x / d))) x

/-- The window below zero: `e^{x/c}` on `x < 0`. -/
noncomputable def winLo (c x : ℝ) : ℝ≥0∞ :=
  (Set.Iio 0).indicator (fun x => ENNReal.ofReal (exp (x / c))) x

theorem measurable_winHi (d : ℝ) : Measurable (winHi d) :=
  (by fun_prop : Measurable fun x : ℝ => ENNReal.ofReal (exp (-x / d))).indicator measurableSet_Ici

theorem measurable_winLo (c : ℝ) : Measurable (winLo c) :=
  (by fun_prop : Measurable fun x : ℝ => ENNReal.ofReal (exp (x / c))).indicator measurableSet_Iio

/-- `A = F − u`: with `X_k` high, the statistic loses the window above zero. -/
theorem sumLaw_hi_add_winHi {d : ℝ} (hd : 0 < d) (m : Multiset ℝ) :
    sumLaw ((-d) ::ₘ m) (Set.Ici 0) + ∫⁻ x, winHi d x ∂(sumLaw m) = sumLaw m (Set.Ici 0) := by
  rw [sumLaw_cons_Ici, ← lintegral_add_left (measurable_expLaw_Ici_neg _),
    ← lintegral_indicator_one measurableSet_Ici]
  refine lintegral_congr fun x => ?_
  rcases le_or_gt 0 x with hx | hx
  · rw [winHi, Set.indicator_of_mem (show x ∈ Set.Ici 0 from hx),
      Set.indicator_of_mem (show x ∈ Set.Ici 0 from hx), Pi.one_apply,
      ← expLaw_neg_Ici_add hd (neg_nonpos.2 hx)]
  · rw [winHi, Set.indicator_of_notMem (show x ∉ Set.Ici (0 : ℝ) from fun h => (not_le.2 hx) h),
      Set.indicator_of_notMem (show x ∉ Set.Ici (0 : ℝ) from fun h => (not_le.2 hx) h),
      expLaw_neg_Ici_of_nonneg hd (by linarith), add_zero]

/-- `B = F + v`: with `X_k` low, the statistic gains the window below zero. -/
theorem sumLaw_lo_eq {c : ℝ} (hc : 0 < c) (m : Multiset ℝ) :
    sumLaw (c ::ₘ m) (Set.Ici 0) = sumLaw m (Set.Ici 0) + ∫⁻ x, winLo c x ∂(sumLaw m) := by
  rw [sumLaw_cons_Ici, ← lintegral_indicator_one measurableSet_Ici,
    ← lintegral_add_left (measurable_one.indicator measurableSet_Ici)]
  refine lintegral_congr fun x => ?_
  rcases le_or_gt 0 x with hx | hx
  · rw [winLo, Set.indicator_of_mem (show x ∈ Set.Ici 0 from hx), Set.indicator_of_notMem
      (show x ∉ Set.Iio (0 : ℝ) from not_lt.2 hx), add_zero, Pi.one_apply,
      expLaw_Ici_of_pos_of_nonpos hc (neg_nonpos.2 hx)]
  · rw [winLo, Set.indicator_of_notMem (show x ∉ Set.Ici (0 : ℝ) from fun h => (not_le.2 hx) h),
      Set.indicator_of_mem (show x ∈ Set.Iio (0 : ℝ) from hx), zero_add,
      expLaw_Ici_of_pos hc (by linarith)]
    congr 2
    ring

/-- `μ` above `ν` in the likelihood-ratio order, tested on ordered pairs of nonnegative functions:
pairs with `ψ₁(t) ψ₂(s) ≤ ψ₁(s) ψ₂(t)` whenever `t ≤ s`. -/
def LRAbove (μ ν : Measure ℝ) : Prop :=
  ∀ ψ₁ ψ₂ : ℝ → ℝ≥0∞, Measurable ψ₁ → Measurable ψ₂ →
    (∀ s t, t ≤ s → ψ₁ t * ψ₂ s ≤ ψ₁ s * ψ₂ t) →
      (∫⁻ x, ψ₁ x ∂ν) * (∫⁻ x, ψ₂ x ∂μ) ≤ (∫⁻ x, ψ₁ x ∂μ) * (∫⁻ x, ψ₂ x ∂ν)

theorem lintegral_mul_lintegral {μ ν : Measure ℝ} [SFinite ν] {f g : ℝ → ℝ≥0∞} (hf : Measurable f)
    (hg : Measurable g) : (∫⁻ x, f x ∂μ) * (∫⁻ y, g y ∂ν) = ∫⁻ x, ∫⁻ y, f x * g y ∂ν ∂μ := by
  rw [← lintegral_mul_const _ hf]
  exact lintegral_congr fun x => (lintegral_const_mul _ hg).symm

/-- Unit exponential integrals compared through their integrands, for `σ, σ' ≥ 0`. -/
theorem lintegral_mul_le_unitExp {f g h k : ℝ → ℝ≥0∞} (hf : Measurable f) (hg : Measurable g)
    (hh : Measurable h) (hk : Measurable k)
    (hpt : ∀ σ σ', 0 ≤ σ → 0 ≤ σ' → f σ * g σ' ≤ h σ' * k σ) :
    (∫⁻ σ, f σ ∂unitExp) * (∫⁻ σ', g σ' ∂unitExp) ≤
      (∫⁻ σ', h σ' ∂unitExp) * (∫⁻ σ, k σ ∂unitExp) := by
  rw [lintegral_mul_lintegral hf hg, mul_comm, lintegral_mul_lintegral hk hh]
  refine lintegral_mono_ae (ae_unitExp_nonneg.mono fun σ hσ => ?_)
  refine lintegral_mono_ae (ae_unitExp_nonneg.mono fun σ' hσ' => ?_)
  rw [mul_comm (k σ)]
  exact hpt σ σ' hσ hσ'

/-- **The base.** `a E` lies above `−b E` in the likelihood-ratio order: their supports are
ordered. -/
theorem lrAbove_expLaw {a b : ℝ} (ha : 0 ≤ a) (hb : 0 ≤ b) : LRAbove (expLaw a) (expLaw (-b)) := by
  intro ψ₁ ψ₂ h₁ h₂ htp
  rw [lintegral_expLaw _ h₁, lintegral_expLaw _ h₂, lintegral_expLaw _ h₁, lintegral_expLaw _ h₂]
  refine lintegral_mul_le_unitExp (h₁.comp (measurable_const_mul _))
    (h₂.comp (measurable_const_mul _)) (h₁.comp (measurable_const_mul _))
    (h₂.comp (measurable_const_mul _)) fun σ σ' hσ hσ' => ?_
  exact htp _ _ (by nlinarith)

/-- Splitting a unit exponential integral at `δ ≥ 0`, by memorylessness. -/
theorem lintegral_unitExp_split {ψ : ℝ → ℝ≥0∞} (hψ : Measurable ψ) (x r : ℝ) {δ : ℝ} (hδ : 0 ≤ δ) :
    ∫⁻ σ, ψ (x + r * σ) ∂unitExp =
      ∫⁻ σ, (Set.Iio δ).indicator (fun σ => ψ (x + r * σ)) σ ∂unitExp +
        ENNReal.ofReal (exp (-δ)) * ∫⁻ σ, ψ ((x + r * δ) + r * σ) ∂unitExp := by
  have hm : Measurable fun σ => ψ (x + r * σ) := hψ.comp (by fun_prop)
  have hre : (fun σ => ψ ((x + r * δ) + r * σ)) = fun σ => ψ (x + r * (σ + δ)) :=
    funext fun σ => congrArg ψ (by ring)
  rw [hre, ← lintegral_unitExp_Ici hm hδ, ← lintegral_add_left (hm.indicator measurableSet_Iio)]
  refine lintegral_congr fun σ => ?_
  rcases lt_or_ge σ δ with h | h
  · rw [Set.indicator_of_mem (show σ ∈ Set.Iio δ from h), Set.indicator_of_notMem
      (show σ ∉ Set.Ici δ from fun h' => (not_le.2 h) h'), add_zero]
  · rw [Set.indicator_of_notMem (show σ ∉ Set.Iio δ from not_lt.2 h), Set.indicator_of_mem
      (show σ ∈ Set.Ici δ from h), zero_add]

/-- **Adding an exponential keeps a pair ordered.** Smoothing an ordered pair `(ψ₁, ψ₂)` by an
independent `r E` gives an ordered pair: split each integral where the two arguments meet
(memorylessness), and the cross terms compare pointwise. -/
theorem smooth_ordered {ψ₁ ψ₂ : ℝ → ℝ≥0∞} (h₁ : Measurable ψ₁) (h₂ : Measurable ψ₂)
    (htp : ∀ s t, t ≤ s → ψ₁ t * ψ₂ s ≤ ψ₁ s * ψ₂ t) (r : ℝ) (s t : ℝ) (hts : t ≤ s) :
    (∫⁻ σ, ψ₁ (t + r * σ) ∂unitExp) * (∫⁻ σ, ψ₂ (s + r * σ) ∂unitExp) ≤
      (∫⁻ σ, ψ₁ (s + r * σ) ∂unitExp) * (∫⁻ σ, ψ₂ (t + r * σ) ∂unitExp) := by
  have hm₁t : Measurable fun σ => ψ₁ (t + r * σ) := h₁.comp (by fun_prop)
  have hm₂t : Measurable fun σ => ψ₂ (t + r * σ) := h₂.comp (by fun_prop)
  have hm₁s : Measurable fun σ => ψ₁ (s + r * σ) := h₁.comp (by fun_prop)
  have hm₂s : Measurable fun σ => ψ₂ (s + r * σ) := h₂.comp (by fun_prop)
  rcases lt_trichotomy r 0 with hr | rfl | hr
  · -- `r < 0`: split the integrals at `s`, which reach `t` after `δ`
    set δ := (s - t) / (-r) with hδ
    have hδ0 : 0 ≤ δ := div_nonneg (by linarith) (by linarith)
    have hsδ : s + r * δ = t := by rw [hδ]; field_simp [hr.ne]; ring
    have e₁ := lintegral_unitExp_split h₁ s r hδ0
    have e₂ := lintegral_unitExp_split h₂ s r hδ0
    rw [hsδ] at e₁ e₂
    have hcross : (∫⁻ σ, (Set.Iio δ).indicator (fun σ => ψ₂ (s + r * σ)) σ ∂unitExp) *
        (∫⁻ σ, ψ₁ (t + r * σ) ∂unitExp) ≤
        (∫⁻ σ, ψ₂ (t + r * σ) ∂unitExp) *
          (∫⁻ σ, (Set.Iio δ).indicator (fun σ => ψ₁ (s + r * σ)) σ ∂unitExp) := by
      refine lintegral_mul_le_unitExp (hm₂s.indicator measurableSet_Iio) hm₁t hm₂t
        (hm₁s.indicator measurableSet_Iio) fun σ σ' hσ hσ' => ?_
      by_cases hσδ : σ < δ
      · rw [Set.indicator_of_mem (show σ ∈ Set.Iio δ from hσδ), Set.indicator_of_mem
          (show σ ∈ Set.Iio δ from hσδ)]
        have : r * δ < r * σ := mul_lt_mul_of_neg_left hσδ hr
        calc ψ₂ (s + r * σ) * ψ₁ (t + r * σ') = ψ₁ (t + r * σ') * ψ₂ (s + r * σ) := mul_comm _ _
          _ ≤ ψ₁ (s + r * σ) * ψ₂ (t + r * σ') := htp _ _ (by nlinarith)
          _ = ψ₂ (t + r * σ') * ψ₁ (s + r * σ) := mul_comm _ _
      · rw [Set.indicator_of_notMem (show σ ∉ Set.Iio δ from hσδ), zero_mul]
        exact zero_le
    rw [e₁, e₂]
    calc (∫⁻ σ, ψ₁ (t + r * σ) ∂unitExp) *
          ((∫⁻ σ, (Set.Iio δ).indicator (fun σ => ψ₂ (s + r * σ)) σ ∂unitExp) +
            ENNReal.ofReal (exp (-δ)) * ∫⁻ σ, ψ₂ (t + r * σ) ∂unitExp)
        = (∫⁻ σ, (Set.Iio δ).indicator (fun σ => ψ₂ (s + r * σ)) σ ∂unitExp) *
            (∫⁻ σ, ψ₁ (t + r * σ) ∂unitExp) + ENNReal.ofReal (exp (-δ)) *
              ((∫⁻ σ, ψ₁ (t + r * σ) ∂unitExp) * ∫⁻ σ, ψ₂ (t + r * σ) ∂unitExp) := by ring
      _ ≤ (∫⁻ σ, ψ₂ (t + r * σ) ∂unitExp) *
            (∫⁻ σ, (Set.Iio δ).indicator (fun σ => ψ₁ (s + r * σ)) σ ∂unitExp) +
            ENNReal.ofReal (exp (-δ)) *
              ((∫⁻ σ, ψ₁ (t + r * σ) ∂unitExp) * ∫⁻ σ, ψ₂ (t + r * σ) ∂unitExp) :=
          add_le_add hcross le_rfl
      _ = _ := by ring
  · simp only [zero_mul, add_zero, lintegral_const, measure_univ, mul_one]
    exact htp s t hts
  · -- `r > 0`: split the integrals at `t`, which reach `s` after `δ`
    set δ := (s - t) / r with hδ
    have hδ0 : 0 ≤ δ := div_nonneg (by linarith) hr.le
    have htδ : t + r * δ = s := by rw [hδ]; field_simp [hr.ne']; ring
    have e₁ := lintegral_unitExp_split h₁ t r hδ0
    have e₂ := lintegral_unitExp_split h₂ t r hδ0
    rw [htδ] at e₁ e₂
    have hcross : (∫⁻ σ, (Set.Iio δ).indicator (fun σ => ψ₁ (t + r * σ)) σ ∂unitExp) *
        (∫⁻ σ, ψ₂ (s + r * σ) ∂unitExp) ≤
        (∫⁻ σ, ψ₁ (s + r * σ) ∂unitExp) *
          (∫⁻ σ, (Set.Iio δ).indicator (fun σ => ψ₂ (t + r * σ)) σ ∂unitExp) := by
      refine lintegral_mul_le_unitExp (hm₁t.indicator measurableSet_Iio) hm₂s hm₁s
        (hm₂t.indicator measurableSet_Iio) fun σ σ' hσ hσ' => ?_
      by_cases hσδ : σ < δ
      · rw [Set.indicator_of_mem (show σ ∈ Set.Iio δ from hσδ), Set.indicator_of_mem
          (show σ ∈ Set.Iio δ from hσδ)]
        have : r * σ < r * δ := mul_lt_mul_of_pos_left hσδ hr
        exact htp _ _ (by nlinarith)
      · rw [Set.indicator_of_notMem (show σ ∉ Set.Iio δ from hσδ), zero_mul]
        exact zero_le
    rw [e₁, e₂]
    calc ((∫⁻ σ, (Set.Iio δ).indicator (fun σ => ψ₁ (t + r * σ)) σ ∂unitExp) +
            ENNReal.ofReal (exp (-δ)) * ∫⁻ σ, ψ₁ (s + r * σ) ∂unitExp) *
          (∫⁻ σ, ψ₂ (s + r * σ) ∂unitExp)
        = (∫⁻ σ, (Set.Iio δ).indicator (fun σ => ψ₁ (t + r * σ)) σ ∂unitExp) *
            (∫⁻ σ, ψ₂ (s + r * σ) ∂unitExp) + ENNReal.ofReal (exp (-δ)) *
              ((∫⁻ σ, ψ₁ (s + r * σ) ∂unitExp) * ∫⁻ σ, ψ₂ (s + r * σ) ∂unitExp) := by ring
      _ ≤ (∫⁻ σ, ψ₁ (s + r * σ) ∂unitExp) *
            (∫⁻ σ, (Set.Iio δ).indicator (fun σ => ψ₂ (t + r * σ)) σ ∂unitExp) +
            ENNReal.ofReal (exp (-δ)) *
              ((∫⁻ σ, ψ₁ (s + r * σ) ∂unitExp) * ∫⁻ σ, ψ₂ (s + r * σ) ∂unitExp) :=
          add_le_add hcross le_rfl
      _ = _ := by ring

/-- **The step.** Adding an independent `r E` to both keeps the likelihood-ratio order. -/
theorem LRAbove.conv {μ ν : Measure ℝ} [SFinite μ] [SFinite ν] (h : LRAbove μ ν) (r : ℝ) :
    LRAbove (μ ∗ expLaw r) (ν ∗ expLaw r) := by
  intro ψ₁ ψ₂ h₁ h₂ htp
  have hsm : ∀ {ψ : ℝ → ℝ≥0∞}, Measurable ψ → Measurable fun x => ∫⁻ σ, ψ (x + r * σ) ∂unitExp :=
    fun hψ =>
      (hψ.comp (by fun_prop : Measurable fun p : ℝ × ℝ => p.1 + r * p.2)).lintegral_prod_right'
  have hconv : ∀ {κ : Measure ℝ} {ψ : ℝ → ℝ≥0∞}, Measurable ψ →
      ∫⁻ x, ψ x ∂(κ ∗ expLaw r) = ∫⁻ x, ∫⁻ σ, ψ (x + r * σ) ∂unitExp ∂κ := by
    intro κ ψ hψ
    rw [Measure.lintegral_conv hψ]
    exact lintegral_congr fun x => lintegral_expLaw r (hψ.comp (measurable_const_add x))
  rw [hconv h₁, hconv h₂, hconv h₁, hconv h₂]
  exact h _ _ (hsm h₁) (hsm h₂) fun s t hts => smooth_ordered h₁ h₂ htp r s t hts

/-- **The likelihood-ratio order along an edge** (Vlassis and Thomas's (22), without densities):
`Y + aE` lies above `Y − bE` for every sum `Y` of independent scaled exponentials, each of
log-concave law. -/
theorem lrAbove_edge {a b : ℝ} (ha : 0 ≤ a) (hb : 0 ≤ b) (Y : Multiset ℝ) :
    LRAbove (sumLaw (a ::ₘ Y)) (sumLaw ((-b) ::ₘ Y)) := by
  induction Y using Multiset.induction_on with
  | empty =>
    rw [sumLaw_cons, sumLaw_cons, sumLaw_zero, Measure.dirac_zero_conv, Measure.dirac_zero_conv]
    exact lrAbove_expLaw ha hb
  | cons r Y ih =>
    rw [Multiset.cons_swap, Multiset.cons_swap (-b), sumLaw_cons r, sumLaw_cons r]
    exact ih.conv r

theorem lintegral_winHi_ne_top {d : ℝ} (hd : 0 < d) (m : Multiset ℝ) :
    ∫⁻ x, winHi d x ∂(sumLaw m) ≠ ⊤ := by
  intro h
  have := sumLaw_hi_add_winHi hd m
  rw [h, add_top] at this
  exact measure_ne_top _ _ this.symm

theorem lintegral_winLo_ne_top {c : ℝ} (hc : 0 < c) (m : Multiset ℝ) :
    ∫⁻ x, winLo c x ∂(sumLaw m) ≠ ⊤ := by
  intro h
  have := sumLaw_lo_eq hc m
  rw [h, add_top] at this
  exact measure_ne_top _ _ this

/-- `F − A`, the mass of the window above zero. -/
theorem stat_sub_hi {d : ℝ} (hd : 0 < d) (m : Multiset ℝ) :
    stat m - stat ((-d) ::ₘ m) = (∫⁻ x, winHi d x ∂(sumLaw m)).toReal := by
  have h := congrArg ENNReal.toReal (sumLaw_hi_add_winHi hd m)
  rw [ENNReal.toReal_add (measure_ne_top _ _) (lintegral_winHi_ne_top hd m)] at h
  simp only [stat, measureReal_def]
  linarith

/-- `B − F`, the mass of the window below zero. -/
theorem stat_lo_sub {c : ℝ} (hc : 0 < c) (m : Multiset ℝ) :
    stat (c ::ₘ m) - stat m = (∫⁻ x, winLo c x ∂(sumLaw m)).toReal := by
  have h := congrArg ENNReal.toReal (sumLaw_lo_eq hc m)
  rw [ENNReal.toReal_add (measure_ne_top _ _) (lintegral_winLo_ne_top hc m)] at h
  simp only [stat, measureReal_def]
  linarith

/-- `w = B − A > 0`: the two windows together are positive everywhere. -/
theorem stat_window_pos {c d : ℝ} (hc : 0 < c) (hd : 0 < d) (m : Multiset ℝ) :
    0 < stat (c ::ₘ m) - stat ((-d) ::ₘ m) := by
  have e : stat (c ::ₘ m) - stat ((-d) ::ₘ m) =
      (∫⁻ x, winHi d x + winLo c x ∂(sumLaw m)).toReal := by
    rw [lintegral_add_left (measurable_winHi d), ENNReal.toReal_add (lintegral_winHi_ne_top hd m)
      (lintegral_winLo_ne_top hc m), ← stat_sub_hi hd, ← stat_lo_sub hc]
    ring
  rw [e]
  refine ENNReal.toReal_pos ?_ (by
    rw [lintegral_add_left (measurable_winHi d)]
    exact ENNReal.add_ne_top.2 ⟨lintegral_winHi_ne_top hd m, lintegral_winLo_ne_top hc m⟩)
  rw [Ne, lintegral_eq_zero_iff (f := fun x => winHi d x + winLo c x)
    ((measurable_winHi d).add (measurable_winLo c))]
  intro h0
  obtain ⟨x, hx⟩ := h0.exists
  have hpos : 0 < winHi d x + winLo c x := by
    rcases le_or_gt 0 x with h | h
    · rw [winHi, Set.indicator_of_mem (show x ∈ Set.Ici 0 from h)]
      exact lt_of_lt_of_le (ENNReal.ofReal_pos.2 (exp_pos _)) le_self_add
    · rw [winLo, Set.indicator_of_mem (show x ∈ Set.Iio 0 from h)]
      exact lt_of_lt_of_le (ENNReal.ofReal_pos.2 (exp_pos _)) le_add_self
  have hx' : winHi d x + winLo c x = 0 := hx
  rw [hx'] at hpos
  exact lt_irrefl _ hpos

/-- **The likelihood-ratio inequality along an edge**, `u₋ v₊ ≤ u₊ v₋`. -/
theorem window_cross {a b c d : ℝ} (ha : 0 ≤ a) (hb : 0 ≤ b) (hc : 0 < c) (hd : 0 < d)
    (Y : Multiset ℝ) :
    (stat ((-b) ::ₘ Y) - stat ((-d) ::ₘ (-b) ::ₘ Y)) * (stat (c ::ₘ a ::ₘ Y) - stat (a ::ₘ Y)) ≤
      (stat (a ::ₘ Y) - stat ((-d) ::ₘ a ::ₘ Y)) *
        (stat (c ::ₘ (-b) ::ₘ Y) - stat ((-b) ::ₘ Y)) := by
  rw [stat_sub_hi hd, stat_sub_hi hd, stat_lo_sub hc, stat_lo_sub hc, ← ENNReal.toReal_mul,
    ← ENNReal.toReal_mul]
  refine ENNReal.toReal_mono (ENNReal.mul_ne_top (lintegral_winHi_ne_top hd _)
    (lintegral_winLo_ne_top hc _)) ?_
  refine lrAbove_edge ha hb Y _ _ (measurable_winHi d) (measurable_winLo c) fun s t hts => ?_
  rcases lt_or_ge t 0 with ht | ht
  · rw [winHi, Set.indicator_of_notMem (show t ∉ Set.Ici (0 : ℝ) from fun h => (not_le.2 ht) h),
      zero_mul]
    exact zero_le
  · rw [winLo, Set.indicator_of_notMem (show s ∉ Set.Iio (0 : ℝ) from not_lt.2 (ht.trans hts)),
      mul_zero]
    exact zero_le

/-- The share of the window above zero in both windows: `θ = (F − A)/(B − A)`. -/
noncomputable def windowShare (c d : ℝ) (m : Multiset ℝ) : ℝ :=
  (stat m - stat ((-d) ::ₘ m)) / (stat (c ::ₘ m) - stat ((-d) ::ₘ m))

/-- **`θ` falls along an edge** (Vlassis and Thomas's (15)): moving one coefficient from `a > 0`
to `−b < 0` lowers the share of the window above zero. -/
theorem windowShare_edge {a b c d : ℝ} (ha : 0 ≤ a) (hb : 0 ≤ b) (hc : 0 < c) (hd : 0 < d)
    (Y : Multiset ℝ) : windowShare c d ((-b) ::ₘ Y) ≤ windowShare c d (a ::ₘ Y) := by
  unfold windowShare
  rw [div_le_div_iff₀ (stat_window_pos hc hd _) (stat_window_pos hc hd _)]
  have h := window_cross ha hb hc hd Y
  nlinarith

/-- **The transfer coefficient's sign** (Vlassis and Thomas's (14) and (31)): along an edge
`Y + aE`, `Y − bE`, with `c ≤ a`, `η = (1 − θ₋)(A₊ − A₋) − p (F₊ − F₋) ≥ 0`, `p = c/(c + d)`.
The race and same-sign identities give the flux identities, and with them
`(c + d) η = (a − c) w₊ (θ₊ − θ₋)`. -/
theorem transfer_nonneg {a b c d : ℝ} (ha : 0 < a) (hb : 0 < b) (hc : 0 < c) (hd : 0 < d)
    (hca : c ≤ a) (Y : Multiset ℝ) :
    0 ≤ (1 - windowShare c d ((-b) ::ₘ Y)) *
        (stat ((-d) ::ₘ a ::ₘ Y) - stat ((-d) ::ₘ (-b) ::ₘ Y)) -
      c / (c + d) * (stat (a ::ₘ Y) - stat ((-b) ::ₘ Y)) := by
  have hR1 := stat_race ha hd Y
  have hS1 := stat_same_neg hb hd Y
  have hR2 := stat_race hc hb Y
  have hS2 := stat_same_pos ha hc Y
  rw [Multiset.cons_swap] at hR1
  rw [Multiset.cons_swap] at hS1
  rw [Multiset.cons_swap] at hS2
  have hθp := windowShare_edge ha.le hb.le hc hd Y
  have hwp := stat_window_pos hc hd (a ::ₘ Y)
  have hwm := stat_window_pos hc hd ((-b) ::ₘ Y)
  set Fp := stat (a ::ₘ Y)
  set Fm := stat ((-b) ::ₘ Y)
  set Ap := stat ((-d) ::ₘ a ::ₘ Y)
  set Am := stat ((-d) ::ₘ (-b) ::ₘ Y)
  set Bp := stat (c ::ₘ a ::ₘ Y)
  set Bm := stat (c ::ₘ (-b) ::ₘ Y)
  set θp := windowShare c d (a ::ₘ Y) with hθpdef
  set θm := windowShare c d ((-b) ::ₘ Y) with hθmdef
  have hp : θp * (Bp - Ap) = Fp - Ap := by
    rw [hθpdef, windowShare]; exact div_mul_cancel₀ _ hwp.ne'
  have hm : θm * (Bm - Am) = Fm - Am := by
    rw [hθmdef, windowShare]; exact div_mul_cancel₀ _ hwm.ne'
  have Φ1 : d * (Ap - Am) = a * (Fp - Ap) + b * (Fm - Am) := by linear_combination hR1 + hS1
  have Φ2 : c * (Bp - Bm) = a * (Bp - Fp) + b * (Bm - Fm) := by linear_combination -hR2 - hS2
  have hkey : ((1 - θm) * (Ap - Am) - c / (c + d) * (Fp - Fm)) * (c + d) =
      (a - c) * (Bp - Ap) * (θp - θm) := by
    have hcd : c / (c + d) * (c + d) = c := div_mul_cancel₀ _ (by linarith)
    linear_combination (1 - θm) * Φ1 + (1 - θm) * Φ2 - (Fp - Fm) * hcd - Φ2 - (a - c) * hp -
      (b + c) * hm
  have hcd0 : 0 < c + d := by linarith
  have := mul_nonneg (mul_nonneg (sub_nonneg.2 hca) hwp.le) (sub_nonneg.2 hθp)
  rw [← hkey] at this
  exact nonneg_of_mul_nonneg_left this hcd0

/-! #### The chain induction (Vlassis and Thomas's §2 and Lemma 6) -/

theorem sum_range_ite_lt_sub (θ : ℕ → ℝ) {j m : ℕ} (hj : j ≤ m) :
    ∑ J ∈ Finset.range (m + 1), (if j < J then θ J - θ (J + 1) else 0) = θ (j + 1) - θ (m + 1) := by
  induction m, hj using Nat.le_induction with
  | base =>
    rw [Finset.sum_eq_zero fun J hJ => if_neg (by simp only [Finset.mem_range] at hJ; omega)]
    ring
  | succ m _ ih =>
    rw [Finset.sum_range_succ, ih, if_pos (by omega)]
    ring

theorem sum_range_ite_le_sub (θ : ℕ → ℝ) {j m : ℕ} (hj : j ≤ m) :
    ∑ J ∈ Finset.range (m + 1), (if J ≤ j then θ J - θ (J + 1) else 0) = θ 0 - θ (j + 1) := by
  induction m, hj using Nat.le_induction with
  | base =>
    rw [Finset.sum_congr rfl fun J hJ => if_pos (by simp only [Finset.mem_range] at hJ; omega),
      Finset.sum_range_sub']
  | succ m _ ih =>
    rw [Finset.sum_range_succ, ih, if_neg (by omega), add_zero]

/-- **The insertion's algebra** (Vlassis and Thomas's (27)–(29)). With
`θ_j (B_j − A_j) = F_j − A_j`, `θ` falling from `1` to `0`, and every transfer `η_j` nonnegative,
the `λ`-average of the inserted chains' means beats the independent reveal's, so one inserted chain
does. -/
theorem insertion_algebra (m : ℕ) (F A B θ gC gH : ℕ → ℝ) (p : ℝ)
    (hw : ∀ j ≤ m + 1, θ j * (B j - A j) = F j - A j)
    (hθ0 : θ 0 = 1) (hθm : θ (m + 1) = 0) (hθ : ∀ j ≤ m, θ (j + 1) ≤ θ j)
    (hη : ∀ j ≤ m, 0 ≤ (1 - θ (j + 1)) * (A j - A (j + 1)) - p * (F j - F (j + 1)))
    (hg : ∀ j ≤ m, gC j ≤ gH j) :
    ∃ J ≤ m, ∑ j ∈ Finset.range (m + 1), (F j - F (j + 1)) * ((1 - p) * gC j + p * gH j) ≤
      ∑ j ∈ Finset.range (m + 1),
        ((if j < J then B j - B (j + 1) else if j = J then B J - A J else 0) * gC j +
          (if J ≤ j then A j - A (j + 1) else 0) * gH j) := by
  set E : ℕ → ℝ := fun J => ∑ j ∈ Finset.range (m + 1),
    ((if j < J then B j - B (j + 1) else if j = J then B J - A J else 0) * gC j +
      (if J ≤ j then A j - A (j + 1) else 0) * gH j) with hE
  obtain ⟨J, hJ, hmax⟩ := Finset.exists_max_image (Finset.range (m + 1)) E ⟨0, by simp⟩
  refine ⟨J, by simp only [Finset.mem_range] at hJ; omega, ?_⟩
  have hlam : ∀ J ∈ Finset.range (m + 1), 0 ≤ θ J - θ (J + 1) := fun J hJ =>
    sub_nonneg.2 (hθ J (by simp only [Finset.mem_range] at hJ; omega))
  have hsumlam : ∑ J ∈ Finset.range (m + 1), (θ J - θ (J + 1)) = 1 := by
    rw [Finset.sum_range_sub', hθ0, hθm, sub_zero]
  -- the λ-average, coefficient by coefficient
  have havg : ∑ J ∈ Finset.range (m + 1), (θ J - θ (J + 1)) * E J =
      ∑ j ∈ Finset.range (m + 1),
        ((θ (j + 1) * (B j - B (j + 1)) + (θ j - θ (j + 1)) * (B j - A j)) * gC j +
          ((1 - θ (j + 1)) * (A j - A (j + 1))) * gH j) := by
    simp only [hE, Finset.mul_sum]
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl fun j hj => ?_
    have hjm : j ≤ m := by simp only [Finset.mem_range] at hj; omega
    have hpt : ∀ J, (θ J - θ (J + 1)) *
        ((if j < J then B j - B (j + 1) else if j = J then B J - A J else 0) * gC j +
          (if J ≤ j then A j - A (j + 1) else 0) * gH j) =
        ((if j < J then θ J - θ (J + 1) else 0) * (B j - B (j + 1)) +
          (if J = j then θ J - θ (J + 1) else 0) * (B j - A j)) * gC j +
        ((if J ≤ j then θ J - θ (J + 1) else 0) * (A j - A (j + 1))) * gH j := by
      intro J
      split_ifs <;> (try subst_vars) <;> first | omega | ring
    rw [Finset.sum_congr rfl fun J _ => hpt J]
    simp only [Finset.sum_add_distrib, ← Finset.sum_mul]
    rw [sum_range_ite_lt_sub θ hjm, sum_range_ite_le_sub θ hjm, Finset.sum_ite_eq', if_pos hj,
      hθm, hθ0]
    ring
  calc ∑ j ∈ Finset.range (m + 1), (F j - F (j + 1)) * ((1 - p) * gC j + p * gH j)
      ≤ ∑ J ∈ Finset.range (m + 1), (θ J - θ (J + 1)) * E J := by
        rw [havg]
        refine Finset.sum_le_sum fun j hj => ?_
        have hjm : j ≤ m := by simp only [Finset.mem_range] at hj; omega
        have key : ((θ (j + 1) * (B j - B (j + 1)) + (θ j - θ (j + 1)) * (B j - A j)) * gC j +
            ((1 - θ (j + 1)) * (A j - A (j + 1))) * gH j) -
            (F j - F (j + 1)) * ((1 - p) * gC j + p * gH j) =
            ((1 - θ (j + 1)) * (A j - A (j + 1)) - p * (F j - F (j + 1))) * (gH j - gC j) := by
          linear_combination (gC j) * (hw j (by omega)) - (gC j) * (hw (j + 1) (by omega))
        have := mul_nonneg (hη j hjm) (sub_nonneg.2 (hg j hjm))
        linarith
    _ ≤ ∑ J' ∈ Finset.range (m + 1), (θ J' - θ (J' + 1)) * E J :=
        Finset.sum_le_sum fun J' hJ' => mul_le_mul_of_nonneg_left (hmax J' hJ') (hlam J' hJ')
    _ = E J := by rw [← Finset.sum_mul, hsumlam, one_mul]

/-- An inserted chain's mean, read by the states before and after the insertion point. -/
theorem chain_insert_sum (A B gC gH : ℕ → ℝ) {J m : ℕ} (hJ : J ≤ m) :
    ∑ k ∈ Finset.range (m + 2),
        ((if k ≤ J then B k else A (k - 1)) - (if k + 1 ≤ J then B (k + 1) else A k)) *
          (if k ≤ J then gC k else gH (k - 1)) =
      ∑ j ∈ Finset.range (m + 1),
        ((if j < J then B j - B (j + 1) else if j = J then B J - A J else 0) * gC j +
          (if J ≤ j then A j - A (j + 1) else 0) * gH j) := by
  induction m, hJ using Nat.le_induction with
  | base =>
    rw [Finset.sum_range_succ, Finset.sum_range_succ, Finset.sum_range_succ (n := J)]
    rw [Finset.sum_congr rfl (g := fun j => (if j < J then B j - B (j + 1) else
        if j = J then B J - A J else 0) * gC j + (if J ≤ j then A j - A (j + 1) else 0) * gH j)
      fun k hk => by
        simp only [Finset.mem_range] at hk
        split_ifs <;> first | omega | ring]
    rw [if_pos le_rfl, if_neg (by omega), if_pos le_rfl, if_neg (lt_irrefl _), if_pos rfl,
      if_pos le_rfl, if_neg (by omega), if_neg (by omega), if_neg (by omega),
      show J + 1 - 1 = J by omega]
    ring
  | succ m hJm ih =>
    rw [Finset.sum_range_succ, ih, Finset.sum_range_succ (n := m + 1), if_neg (by omega),
      if_neg (by omega), if_neg (by omega), if_neg (by omega), if_neg (by omega),
      if_pos (by omega), show m + 2 - 1 = m + 1 by omega]
    ring

variable {ι : Type*} [DecidableEq ι]

/-- A unit's coefficient in the two-point statistic: `−β` when high, `γ` when low. -/
def coef (γ β : ι → ℝ) (S : Finset ι) (i : ι) : ℝ := if i ∈ S then -β i else γ i

/-- The two-point statistic of the units in `P` with high set `S` (Vlassis and Thomas's
`K_m(S)`): `P(E₀ + ∑_{i ∈ P \ S} γ_i E_i − ∑_{i ∈ S} β_i E_i ≥ 0)`. -/
noncomputable def cubeK (γ β : ι → ℝ) (P S : Finset ι) : ℝ :=
  stat (1 ::ₘ P.val.map (coef γ β S))

/-- A maximal chain through `P`: `S 0 = ∅`, `S |P| = P`, one unit of `P` added at each step. -/
def IsChain (P : Finset ι) (S : ℕ → Finset ι) : Prop :=
  S 0 = ∅ ∧ S P.card = P ∧ ∀ j < P.card, ∃ e ∈ P, e ∉ S j ∧ S (j + 1) = insert e (S j)

/-- The mean of `g` under a chain's measure (Vlassis and Thomas's `ν_C`): mass
`K(S_j) − K(S_{j+1})` on `S_j`, and `K(S_m)` on the last state. -/
noncomputable def chainMean (γ β : ι → ℝ) (P : Finset ι) (S : ℕ → Finset ι)
    (g : Finset ι → ℝ) : ℝ :=
  ∑ j ∈ Finset.range (P.card + 1),
    (cubeK γ β P (S j) - if j < P.card then cubeK γ β P (S (j + 1)) else 0) * g (S j)

theorem IsChain.subset {P : Finset ι} {S : ℕ → Finset ι} (hS : IsChain P S) :
    ∀ j ≤ P.card, S j ⊆ P := by
  intro j
  induction j with
  | zero => intro _; rw [hS.1]; exact Finset.empty_subset _
  | succ j ih =>
    intro hj
    obtain ⟨e, he, -, hSe⟩ := hS.2.2 j (by omega)
    rw [hSe]
    exact Finset.insert_subset he (ih (by omega))

/-- `K` falls as the high set grows. -/
theorem cubeK_anti {γ β : ι → ℝ} (P : Finset ι) (hγβ : ∀ i ∈ P, -β i ≤ γ i) {S S' : Finset ι}
    (h : S ⊆ S') : cubeK γ β P S' ≤ cubeK γ β P S := by
  unfold cubeK
  rw [← Multiset.singleton_add, ← Multiset.singleton_add]
  refine stat_map_mono {1} P.val fun i hi => ?_
  unfold coef
  by_cases hS : i ∈ S
  · rw [if_pos (h hS), if_pos hS]
  · by_cases hS' : i ∈ S'
    · rw [if_pos hS', if_neg hS]; exact hγβ i hi
    · rw [if_neg hS', if_neg hS]

theorem coef_insert_of_ne {γ β : ι → ℝ} {S : Finset ι} {e x : ι} (hx : x ≠ e) :
    coef γ β (insert e S) x = coef γ β S x := by
  unfold coef
  simp only [Finset.mem_insert, hx, false_or]

theorem map_coef_insert {γ β : ι → ℝ} {P : Finset ι} {i : ι} (hi : i ∉ P) (T : Finset ι) :
    (insert i P).val.map (coef γ β T) = coef γ β T i ::ₘ P.val.map (coef γ β T) := by
  rw [Finset.insert_val_of_notMem hi, Multiset.map_cons]

theorem map_coef_erase {γ β : ι → ℝ} {P : Finset ι} {e : ι} (he : e ∈ P) (T : Finset ι) :
    P.val.map (coef γ β T) = coef γ β T e ::ₘ (P.erase e).val.map (coef γ β T) := by
  rw [Finset.erase_val, ← Multiset.map_cons, Multiset.cons_erase he]

/-- The independent two-point law on the subsets of `s`: unit `i` high with chance `p i`. -/
def bern (p : ι → ℝ) (s T : Finset ι) : ℝ := (∏ i ∈ T, p i) * ∏ i ∈ s \ T, (1 - p i)

theorem bern_nonneg {p : ι → ℝ} (hp : ∀ i, 0 ≤ p i ∧ p i ≤ 1) (s T : Finset ι) : 0 ≤ bern p s T :=
  mul_nonneg (Finset.prod_nonneg fun i _ => (hp i).1)
    (Finset.prod_nonneg fun i _ => sub_nonneg.2 (hp i).2)

/-- Revealing one more independent unit. -/
theorem sum_bern_insert (p : ι → ℝ) {s : Finset ι} {i : ι} (hi : i ∉ s) (f : Finset ι → ℝ) :
    ∑ T ∈ (insert i s).powerset, bern p (insert i s) T * f T =
      (1 - p i) * ∑ T ∈ s.powerset, bern p s T * f T +
        p i * ∑ T ∈ s.powerset, bern p s T * f (insert i T) := by
  rw [Finset.sum_powerset_insert hi, Finset.mul_sum, Finset.mul_sum]
  congr 1
  · refine Finset.sum_congr rfl fun T hT => ?_
    have hiT : i ∉ T := fun h => hi (Finset.mem_powerset.1 hT h)
    rw [bern, bern, Finset.insert_sdiff_of_notMem _ hiT, Finset.prod_insert
      (fun h => hi (Finset.mem_sdiff.1 h).1)]
    ring
  · refine Finset.sum_congr rfl fun T hT => ?_
    have hiT : i ∉ T := fun h => hi (Finset.mem_powerset.1 hT h)
    rw [bern, bern, Finset.prod_insert hiT, Finset.insert_sdiff_insert,
      Finset.sdiff_insert_of_notMem hi]
    ring

/-- The chance a two-point unit is high: `γ/(γ + β)`, the mean-one weight on `1 + β`. -/
noncomputable def highChance (γ β : ι → ℝ) (i : ι) : ℝ := γ i / (γ i + β i)

omit [DecidableEq ι] in
theorem highChance_mem {γ β : ι → ℝ} (hγ : ∀ e, 0 ≤ γ e) (hβ : ∀ e, 0 < β e) (i : ι) :
    0 ≤ highChance γ β i ∧ highChance γ β i ≤ 1 :=
  ⟨div_nonneg (hγ i) (by linarith [hγ i, hβ i]),
    (div_le_one (by linarith [hγ i, hβ i])).2 (by linarith [hβ i])⟩

/-- **The local insertion** (Vlassis and Thomas's Lemma 6). A chain through `P`, a new unit `i` of
least `γ`, and an increasing `g`: some insertion of `i` into the chain has a mean of `g` at least
that of the chain with `i` revealed independently. With `γ_i = 0` the unit is never high, and it
goes last. -/
theorem local_insertion {γ β : ι → ℝ} {P : Finset ι} {S : ℕ → Finset ι} (hS : IsChain P S)
    {i : ι} (hi : i ∉ P) (hγ : ∀ e, 0 ≤ γ e) (hγ1 : γ i ≤ 1) (hβ : ∀ e, 0 < β e)
    (hmin : ∀ e ∈ P, γ i ≤ γ e) {g : Finset ι → ℝ} (hg : ∀ U V, U ⊆ V → g U ≤ g V) :
    ∃ S', IsChain (insert i P) S' ∧
      ∑ j ∈ Finset.range (P.card + 1),
          (cubeK γ β P (S j) - if j < P.card then cubeK γ β P (S (j + 1)) else 0) *
            ((1 - highChance γ β i) * g (S j) + highChance γ β i * g (insert i (S j))) ≤
        chainMean γ β (insert i P) S' g := by
  classical
  set m := P.card with hm
  set c := γ i with hcdef
  set d := β i with hddef
  have hd : 0 < d := hβ i
  set mult : ℕ → Multiset ℝ := fun j => if j ≤ m then 1 ::ₘ P.val.map (coef γ β (S j))
    else (-1) ::ₘ P.val.map (fun e => -β e) with hmult
  set F : ℕ → ℝ := fun j => stat (mult j) with hFdef
  set A : ℕ → ℝ := fun j => stat ((-d) ::ₘ mult j) with hAdef
  set B : ℕ → ℝ := fun j => stat (c ::ₘ mult j) with hBdef
  have hiS : ∀ j ≤ m, i ∉ S j := fun j hj h => hi (hS.subset j hj h)
  have hmultj : ∀ j ≤ m, mult j = 1 ::ₘ P.val.map (coef γ β (S j)) := fun j hj => by
    simp only [hmult, if_pos hj]
  have hF : ∀ j ≤ m, cubeK γ β P (S j) = F j := fun j hj => by
    simp only [hFdef, hmultj j hj, cubeK]
  have hB : ∀ j ≤ m, cubeK γ β (insert i P) (S j) = B j := fun j hj => by
    simp only [hBdef, hmultj j hj, cubeK, map_coef_insert hi]
    rw [show coef γ β (S j) i = c by simp [coef, hiS j hj, hcdef], Multiset.cons_swap]
  have hA : ∀ j ≤ m, cubeK γ β (insert i P) (insert i (S j)) = A j := fun j hj => by
    simp only [hAdef, hmultj j hj, cubeK, map_coef_insert hi]
    have hmap : P.val.map (coef γ β (insert i (S j))) = P.val.map (coef γ β (S j)) :=
      Multiset.map_congr rfl fun x hx => coef_insert_of_ne (fun h => hi (by rw [← h]; exact hx))
    rw [show coef γ β (insert i (S j)) i = -d by simp [coef, hddef], hmap, Multiset.cons_swap]
  have hmultm : mult (m + 1) = (-1) ::ₘ P.val.map (fun e => -β e) := by
    simp only [hmult, if_neg (show ¬ m + 1 ≤ m by omega)]
  have hnonpos : ∀ x ∈ P.val.map (fun e => -β e), x ≤ 0 := fun x hx => by
    obtain ⟨e, -, rfl⟩ := Multiset.mem_map.1 hx
    linarith [hβ e]
  have hFm : F (m + 1) = 0 := by
    simp only [hFdef, hmultm]
    exact stat_eq_zero_of_nonpos (by norm_num) hnonpos
  have hAm : A (m + 1) = 0 := by
    simp only [hAdef, hmultm]
    refine stat_eq_zero_of_nonpos (neg_lt_zero.2 hd) fun x hx => ?_
    rcases Multiset.mem_cons.1 hx with rfl | hx
    · norm_num
    · exact hnonpos x hx
  -- the inserted chains
  set S' : ℕ → ℕ → Finset ι := fun J k => if k ≤ J then S k else insert i (S (k - 1)) with hS'
  have hcard : (insert i P).card = m + 1 := Finset.card_insert_of_notMem hi
  have hchain : ∀ J ≤ m, IsChain (insert i P) (S' J) := by
    intro J hJ
    refine ⟨by simp [hS', hS.1], ?_, fun k hk => ?_⟩
    · rw [hcard]
      simp only [hS', if_neg (show ¬ m + 1 ≤ J by omega), show m + 1 - 1 = m by omega]
      rw [hS.2.1]
    · rw [hcard] at hk
      rcases lt_trichotomy k J with h | rfl | h
      · obtain ⟨e, he, heS, hSe⟩ := hS.2.2 k (by omega)
        refine ⟨e, Finset.mem_insert_of_mem he, ?_, ?_⟩
        · simp only [hS', if_pos h.le]; exact heS
        · simp only [hS', if_pos h.le, if_pos (show k + 1 ≤ J by omega)]
          exact hSe
      · refine ⟨i, Finset.mem_insert_self _ _, ?_, ?_⟩
        · simp only [hS', if_pos le_rfl]; exact hiS k (by omega)
        · simp only [hS', if_pos le_rfl, if_neg (show ¬ k + 1 ≤ k by omega),
            show k + 1 - 1 = k by omega]
      · obtain ⟨e, he, heS, hSe⟩ := hS.2.2 (k - 1) (by omega)
        refine ⟨e, Finset.mem_insert_of_mem he, ?_, ?_⟩
        · simp only [hS', if_neg (show ¬ k ≤ J by omega)]
          rw [Finset.mem_insert, not_or]
          exact ⟨fun h' => hi (h' ▸ he), heS⟩
        · simp only [hS', if_neg (show ¬ k ≤ J by omega), if_neg (show ¬ k + 1 ≤ J by omega),
            show k + 1 - 1 = k - 1 + 1 by omega, hSe, Finset.insert_comm]
  set gC : ℕ → ℝ := fun j => g (S j) with hgC
  set gH : ℕ → ℝ := fun j => g (insert i (S j)) with hgH
  have hmean : ∀ J ≤ m, chainMean γ β (insert i P) (S' J) g =
      ∑ j ∈ Finset.range (m + 1),
        ((if j < J then B j - B (j + 1) else if j = J then B J - A J else 0) * gC j +
          (if J ≤ j then A j - A (j + 1) else 0) * gH j) := by
    intro J hJ
    rw [← chain_insert_sum A B gC gH hJ, chainMean, hcard]
    refine Finset.sum_congr rfl fun k hk => ?_
    simp only [Finset.mem_range] at hk
    have hK : cubeK γ β (insert i P) (S' J k) = if k ≤ J then B k else A (k - 1) := by
      by_cases hkJ : k ≤ J
      · simp only [hS', if_pos hkJ]; exact hB k (by omega)
      · simp only [hS', if_neg hkJ]; exact hA (k - 1) (by omega)
    have hK' : (if k < m + 1 then cubeK γ β (insert i P) (S' J (k + 1)) else 0) =
        if k + 1 ≤ J then B (k + 1) else A k := by
      by_cases hkm : k < m + 1
      · rw [if_pos hkm]
        by_cases hkJ : k + 1 ≤ J
        · simp only [hS', if_pos hkJ]; exact hB (k + 1) (by omega)
        · simp only [hS', if_neg hkJ, show k + 1 - 1 = k by omega]; exact hA k (by omega)
      · rw [if_neg hkm, if_neg (by omega), show k = m + 1 by omega, hAm]
    have hgk : g (S' J k) = if k ≤ J then gC k else gH (k - 1) := by
      by_cases hkJ : k ≤ J
      · simp only [hS', hgC, if_pos hkJ]
      · simp only [hS', hgH, if_neg hkJ]
    rw [hK, hK', hgk]
  have hhat : ∑ j ∈ Finset.range (m + 1),
      (cubeK γ β P (S j) - if j < m then cubeK γ β P (S (j + 1)) else 0) *
        ((1 - highChance γ β i) * g (S j) + highChance γ β i * g (insert i (S j))) =
      ∑ j ∈ Finset.range (m + 1), (F j - F (j + 1)) *
        ((1 - highChance γ β i) * gC j + highChance γ β i * gH j) := by
    refine Finset.sum_congr rfl fun j hj => ?_
    simp only [Finset.mem_range] at hj
    rw [hF j (by omega)]
    by_cases hjm : j < m
    · rw [if_pos hjm, hF (j + 1) (by omega)]
    · rw [if_neg hjm, show j + 1 = m + 1 by omega, hFm]
  have hgCH : ∀ j ≤ m, gC j ≤ gH j := fun j _ => hg _ _ (Finset.subset_insert _ _)
  rw [hhat]
  by_cases hc0 : c = 0
  · -- `γ_i = 0`: never high; insert it last
    refine ⟨S' m, hchain m le_rfl, ?_⟩
    have hp0 : highChance γ β i = 0 := by simp [highChance, ← hcdef, hc0]
    have hBF : ∀ j, B j = F j := fun j => by simp only [hBdef, hFdef, hc0, stat_cons_zero]
    rw [hmean m le_rfl, hp0, Finset.sum_range_succ, Finset.sum_range_succ]
    refine add_le_add (le_of_eq (Finset.sum_congr rfl fun j hj => ?_)) ?_
    · simp only [Finset.mem_range] at hj
      rw [if_pos hj, if_neg (show ¬ m ≤ j by omega), hBF, hBF]
      ring
    · rw [if_neg (lt_irrefl _), if_pos rfl, if_pos le_rfl, hAm, hFm, hBF]
      have h1 : 0 ≤ A m * (gH m - gC m) := mul_nonneg (stat_nonneg _) (sub_nonneg.2 (hgCH m le_rfl))
      nlinarith [h1]
  · have hc : 0 < c := lt_of_le_of_ne (hγ i) (Ne.symm hc0)
    set θ : ℕ → ℝ := fun j => windowShare c d (mult j) with hθdef
    have hw : ∀ j ≤ m + 1, θ j * (B j - A j) = F j - A j := fun j _ => by
      simp only [hθdef, windowShare, hFdef, hAdef, hBdef]
      exact div_mul_cancel₀ _ (stat_window_pos hc hd _).ne'
    have hθ0 : θ 0 = 1 := by
      have hm0 : mult 0 = 1 ::ₘ P.val.map γ := by
        rw [hmultj 0 (Nat.zero_le _), hS.1]
        exact congrArg _ (Multiset.map_congr rfl fun x _ => by simp [coef])
      have hnn : ∀ r ∈ mult 0, 0 ≤ r := by
        rw [hm0]
        intro r hr
        rcases Multiset.mem_cons.1 hr with rfl | hr
        · norm_num
        · obtain ⟨e, -, rfl⟩ := Multiset.mem_map.1 hr
          exact hγ e
      have hF0 : stat (mult 0) = 1 := stat_eq_one_of_nonneg hnn
      have hB0 : stat (c ::ₘ mult 0) = 1 := stat_eq_one_of_nonneg fun r hr => by
        rcases Multiset.mem_cons.1 hr with rfl | hr
        · exact hc.le
        · exact hnn r hr
      have hw0 := stat_window_pos hc hd (mult 0)
      simp only [hθdef, windowShare]
      rw [hF0, hB0]
      rw [hB0] at hw0
      exact div_self hw0.ne'
    have hθm : θ (m + 1) = 0 := by
      have h1 : stat (mult (m + 1)) = 0 := hFm
      have h2 : stat ((-d) ::ₘ mult (m + 1)) = 0 := hAm
      simp only [hθdef, windowShare, h1, h2, sub_zero, zero_div]
    -- every edge of the chain, and the sentinel's
    have hedge : ∀ j ≤ m, ∃ a b : ℝ, ∃ Y : Multiset ℝ, 0 < a ∧ 0 < b ∧ c ≤ a ∧
        mult j = a ::ₘ Y ∧ mult (j + 1) = (-b) ::ₘ Y := by
      intro j hj
      rcases lt_or_eq_of_le hj with hjm | rfl
      · obtain ⟨e, he, heS, hSe⟩ := hS.2.2 j hjm
        refine ⟨γ e, β e, 1 ::ₘ (P.erase e).val.map (coef γ β (S j)),
          lt_of_lt_of_le hc (hmin e he), hβ e, hmin e he, ?_, ?_⟩
        · rw [hmultj j hj, map_coef_erase he, show coef γ β (S j) e = γ e by simp [coef, heS],
            Multiset.cons_swap]
        · rw [hmultj (j + 1) (by omega), hSe, map_coef_erase he,
            show coef γ β (insert e (S j)) e = -β e by simp [coef],
            Multiset.map_congr rfl fun x hx => coef_insert_of_ne
              (fun h => Finset.notMem_erase e P (by subst h; exact hx)),
            Multiset.cons_swap]
      · refine ⟨1, 1, P.val.map (fun e => -β e), one_pos, one_pos, hγ1, ?_, hmultm⟩
        rw [hmultj _ le_rfl, hS.2.1]
        congr 1
        exact Multiset.map_congr rfl fun x hx => by simp [coef, Finset.mem_val.mp hx]
    have hθ : ∀ j ≤ m, θ (j + 1) ≤ θ j := fun j hj => by
      obtain ⟨a, b, Y, ha, hb, -, h1, h2⟩ := hedge j hj
      simp only [hθdef, h1, h2]
      exact windowShare_edge ha.le hb.le hc hd Y
    have hη : ∀ j ≤ m, 0 ≤ (1 - θ (j + 1)) * (A j - A (j + 1)) -
        highChance γ β i * (F j - F (j + 1)) := fun j hj => by
      obtain ⟨a, b, Y, ha, hb, hca, h1, h2⟩ := hedge j hj
      simp only [hθdef, hAdef, hFdef, h1, h2, highChance, ← hcdef, ← hddef]
      exact transfer_nonneg ha hb hc hd hca Y
    obtain ⟨J, hJ, hle⟩ := insertion_algebra m F A B θ gC gH (highChance γ β i) hw hθ0 hθm hθ hη
      hgCH
    exact ⟨S' J, hchain J hJ, by rw [hmean J hJ]; exact hle⟩

/-- The terminal telescoping: the chain's masses sum to `K(S₀)`. -/
theorem sum_range_chain_mass (K : ℕ → ℝ) (m : ℕ) :
    ∑ j ∈ Finset.range (m + 1), (K j - if j < m then K (j + 1) else 0) = K 0 := by
  induction m generalizing K with
  | zero => simp
  | succ m ih =>
    rw [Finset.sum_range_succ', if_pos (Nat.succ_pos m)]
    have := ih (fun j => K (j + 1))
    simp only [show ∀ j, (j + 1 < m + 1) = (j < m) from fun j => propext Nat.succ_lt_succ_iff]
      at this ⊢
    rw [this]
    ring

/-- **A chain measure puts at most `α` on `{K ≤ α}`** (the end of Vlassis and Thomas's
Proposition 3): the rejected states are a terminal segment, whose masses telescope to `K` at its
first state. -/
theorem chain_terminal {α : ℝ} (m : ℕ) {K : ℕ → ℝ} (hK : ∀ j < m, K (j + 1) ≤ K j) :
    ∑ j ∈ Finset.range (m + 1),
      (K j - if j < m then K (j + 1) else 0) * (if K j ≤ α then 1 else 0) ≤ max α 0 := by
  induction m generalizing K with
  | zero => simp only [zero_add, Finset.range_one, Finset.sum_singleton, lt_irrefl, if_false,
      sub_zero]; split_ifs with h <;> simp [h]
  | succ m ih =>
    by_cases h0 : K 0 ≤ α
    · have hle : ∀ j ≤ m + 1, K j ≤ K 0 := by
        intro j
        induction j with
        | zero => intro _; exact le_rfl
        | succ j ihj => intro hj; exact (hK j (by omega)).trans (ihj (by omega))
      rw [Finset.sum_congr rfl fun j hj => by
        rw [if_pos ((hle j (by simp only [Finset.mem_range] at hj; omega)).trans h0), mul_one],
        sum_range_chain_mass]
      exact h0.trans (le_max_left _ _)
    · rw [Finset.sum_range_succ', if_neg h0, mul_zero, add_zero]
      have := ih (K := fun j => K (j + 1)) fun j hj => hK (j + 1) (by omega)
      simp only [show ∀ j, (j + 1 < m + 1) = (j < m) from fun j => propext Nat.succ_lt_succ_iff]
        at this ⊢
      exact this

/-- **The chain induction** (Vlassis and Thomas's §2): placing the units of `P` one at a time,
least `γ` last, some chain through `P`, with the units of `A` outside `P` still independent, has
an increasing payoff's mean at least its mean under the independent law. -/
theorem chain_dominates {γ β : ι → ℝ} (hγ : ∀ e, 0 ≤ γ e) (hγ1 : ∀ e, γ e ≤ 1)
    (hβ : ∀ e, 0 < β e) (A : Finset ι) {h : Finset ι → ℝ} (hh : ∀ U V, U ⊆ V → h U ≤ h V) :
    ∀ P ⊆ A, ∃ S, IsChain P S ∧
      ∑ T ∈ A.powerset, bern (highChance γ β) A T * h T ≤
        chainMean γ β P S
          (fun U => ∑ T ∈ (A \ P).powerset, bern (highChance γ β) (A \ P) T * h (U ∪ T)) := by
  classical
  have hp := highChance_mem hγ hβ
  intro P
  induction P using Finset.strongInduction with
  | H P ih =>
    intro hPA
    rcases P.eq_empty_or_nonempty with rfl | hne
    · refine ⟨fun _ => ∅, ⟨rfl, rfl, fun j hj => by simp at hj⟩, le_of_eq ?_⟩
      simp only [chainMean, Finset.card_empty, zero_add, Finset.range_one, Finset.sum_singleton,
        lt_irrefl, if_false, sub_zero, Finset.sdiff_empty, Finset.empty_union]
      rw [show cubeK γ β ∅ ∅ = 1 from stat_eq_one_of_nonneg (by simp), one_mul]
    · obtain ⟨i, hiP, hmin⟩ := Finset.exists_min_image P γ hne
      have hiA : i ∈ A := hPA hiP
      obtain ⟨S, hS, hdom⟩ := ih (P.erase i) (Finset.erase_ssubset hiP)
        ((Finset.erase_subset i P).trans hPA)
      set g : Finset ι → ℝ := fun U =>
        ∑ T ∈ (A \ P).powerset, bern (highChance γ β) (A \ P) T * h (U ∪ T) with hgdef
      have hg : ∀ U V, U ⊆ V → g U ≤ g V := fun U V hUV =>
        Finset.sum_le_sum fun T _ => mul_le_mul_of_nonneg_left
          (hh _ _ (Finset.union_subset_union hUV le_rfl)) (bern_nonneg hp _ _)
      have hsd : A \ P.erase i = insert i (A \ P) := Finset.sdiff_erase hiA
      have hiAP : i ∉ A \ P := fun h => (Finset.mem_sdiff.1 h).2 hiP
      have hsplit : ∀ U, ∑ T ∈ (A \ P.erase i).powerset,
          bern (highChance γ β) (A \ P.erase i) T * h (U ∪ T) =
          (1 - highChance γ β i) * g U + highChance γ β i * g (insert i U) := fun U => by
        rw [hsd, sum_bern_insert _ hiAP]
        congr 2
        refine Finset.sum_congr rfl fun T _ => ?_
        rw [Finset.union_insert, Finset.insert_union]
      obtain ⟨S', hS', hle⟩ := local_insertion hS (Finset.notMem_erase i P) hγ (hγ1 i) hβ
        (fun e he => hmin e (Finset.mem_of_mem_erase he)) hg
      rw [Finset.insert_erase hiP] at hS' hle
      refine ⟨S', hS', hdom.trans (le_of_eq_of_le ?_ hle)⟩
      simp only [chainMean]
      exact Finset.sum_congr rfl fun j _ => by rw [hsplit]

/-- **Two-point validity** (Vlassis and Thomas's Proposition 3). Independent mean-one two-point
units, low `1 − γ_i` with `0 ≤ γ_i ≤ 1` and high `1 + β_i` with `β_i > 0`: their statistic is at
most `α` with chance at most `α`. -/
theorem twoPoint_valid [Fintype ι] {γ β : ι → ℝ} (hγ : ∀ e, 0 ≤ γ e) (hγ1 : ∀ e, γ e ≤ 1)
    (hβ : ∀ e, 0 < β e) (α : ℝ) :
    ∑ T ∈ (Finset.univ : Finset ι).powerset, bern (highChance γ β) Finset.univ T *
      (if cubeK γ β Finset.univ T ≤ α then 1 else 0) ≤ max α 0 := by
  classical
  have hγβ : ∀ e ∈ (Finset.univ : Finset ι), -β e ≤ γ e := fun e _ => by linarith [hγ e, hβ e]
  have hh : ∀ U V : Finset ι, U ⊆ V → (if cubeK γ β Finset.univ U ≤ α then (1 : ℝ) else 0) ≤
      (if cubeK γ β Finset.univ V ≤ α then 1 else 0) := fun U V hUV => by
    split_ifs with h1 h2
    · exact le_rfl
    · exact absurd ((cubeK_anti _ hγβ hUV).trans h1) h2
    · norm_num
    · exact le_rfl
  obtain ⟨S, hS, hdom⟩ := chain_dominates hγ hγ1 hβ Finset.univ hh Finset.univ le_rfl
  refine hdom.trans (le_of_eq_of_le ?_ (chain_terminal (K := fun j =>
    cubeK γ β Finset.univ (S j)) (Finset.univ : Finset ι).card fun j hj => ?_))
  · simp only [chainMean, Finset.sdiff_self, Finset.powerset_empty, Finset.sum_singleton, bern,
      Finset.prod_empty, Finset.empty_sdiff, one_mul, Finset.union_empty]
  · obtain ⟨e, -, -, hSe⟩ := hS.2.2 j hj
    exact cubeK_anti _ hγβ (hSe ▸ Finset.subset_insert e (S j))

/-! #### Every finite law, and Gaffke's bound -/

omit [DecidableEq ι] in
theorem sumLaw_map_clm (L : ℝ →L[ℝ] ℝ) (m : Multiset ℝ) :
    sumLaw (m.map L) = (sumLaw m).map L := by
  induction m using Multiset.induction_on with
  | empty => rw [Multiset.map_zero, sumLaw_zero, Measure.map_dirac' L.measurable, map_zero]
  | cons r m ih =>
    rw [Multiset.map_cons, sumLaw_cons, sumLaw_cons, ih, Measure.map_conv_continuousLinearMap]
    congr 1
    rw [expLaw, expLaw, Measure.map_map L.measurable (measurable_const_mul r)]
    congr 1
    funext s
    change L r * s = L (r * s)
    rw [show r * s = s • r by rw [smul_eq_mul, mul_comm], L.map_smul, smul_eq_mul, mul_comm]

omit [DecidableEq ι] in
/-- `stat` is blind to a positive scale. -/
theorem stat_map_mul {c : ℝ} (hc : 0 < c) (m : Multiset ℝ) :
    stat (m.map fun r => c * r) = stat m := by
  have h := sumLaw_map_clm (c • ContinuousLinearMap.id ℝ ℝ) m
  have hL : ((c • ContinuousLinearMap.id ℝ ℝ : ℝ →L[ℝ] ℝ) : ℝ → ℝ) = fun r => c * r := by
    funext r; simp
  rw [hL] at h
  unfold stat
  rw [h, measureReal_def, measureReal_def,
    Measure.map_apply (measurable_const_mul c) measurableSet_Ici]
  congr 2
  ext r
  simp [mul_nonneg_iff_of_pos_left hc]

omit [DecidableEq ι] in
theorem sumLaw_singleton_zero {m : Multiset ℝ} {r : ℝ} (hr : r ∈ m) (hr0 : r ≠ 0) :
    sumLaw m {0} = 0 := by
  obtain ⟨m', rfl⟩ := Multiset.exists_cons_of_mem hr
  rw [sumLaw_cons, ← lintegral_indicator_one (measurableSet_singleton 0),
    Measure.lintegral_conv (measurable_one.indicator (measurableSet_singleton 0))]
  have h : ∀ x : ℝ, ∫⁻ y, ({0} : Set ℝ).indicator 1 (x + y) ∂(expLaw r) = 0 := fun x => by
    rw [show (fun y => ({0} : Set ℝ).indicator 1 (x + y)) = ({-x} : Set ℝ).indicator 1 by
      funext y
      simp only [Set.indicator, Set.mem_singleton_iff, Pi.one_apply]
      congr 1
      exact propext ⟨fun h => by linarith, fun h => by linarith⟩,
      lintegral_indicator_one (measurableSet_singleton _), expLaw_singleton hr0]
  exact (lintegral_congr fun x => h x).trans lintegral_zero

omit [DecidableEq ι] in
/-- `P(Z ≥ 0) + P(Z ≤ 0) = 1` when some coefficient is nonzero: no atom at zero. -/
theorem stat_add_stat_neg {m : Multiset ℝ} {r : ℝ} (hr : r ∈ m) (hr0 : r ≠ 0) :
    stat m + stat (m.map Neg.neg) = 1 := by
  have h := sumLaw_map_clm (-ContinuousLinearMap.id ℝ ℝ) m
  have hL : ((-ContinuousLinearMap.id ℝ ℝ : ℝ →L[ℝ] ℝ) : ℝ → ℝ) = Neg.neg := by
    funext r; simp
  rw [hL] at h
  unfold stat
  rw [h, measureReal_def, measureReal_def, Measure.map_apply measurable_neg measurableSet_Ici,
    show Neg.neg ⁻¹' Set.Ici (0 : ℝ) = Set.Iic 0 by ext x; simp,
    ← ENNReal.toReal_add (measure_ne_top _ _) (measure_ne_top _ _)]
  have hu := measure_union_add_inter (μ := sumLaw m) (Set.Ici (0 : ℝ))
    (measurableSet_Iic (a := (0 : ℝ)))
  rw [Set.Ici_inter_Iic, Set.Icc_self, sumLaw_singleton_zero hr hr0, add_zero,
    show Set.Ici (0 : ℝ) ∪ Set.Iic 0 = Set.univ by ext x; simp, measure_univ] at hu
  rw [← hu, ENNReal.toReal_one]

/-- **Vlassis and Thomas's statistic**, `K(x) = P(∑ x_i D_i ≤ 1)` with `D ~ Dir(1, …, 1)` over
the `x_i` and an added zero; as exponentials, `P(E₀ + ∑ (1 − x_i) E_i ≥ 0)`. -/
noncomputable def gaffkeK {ι : Type*} [Fintype ι] (x : ι → ℝ) : ℝ :=
  stat (1 ::ₘ Finset.univ.val.map fun i => 1 - x i)

theorem gaffkeK_anti {ι : Type*} [Fintype ι] {x y : ι → ℝ} (h : ∀ i, x i ≤ y i) :
    gaffkeK y ≤ gaffkeK x := by
  unfold gaffkeK
  rw [← Multiset.singleton_add, ← Multiset.singleton_add]
  exact stat_map_mono {1} _ fun i _ => by linarith [h i]

/-- **The cube meets `K`**: the two-point outcome with high set `T` has statistic
`cubeK γ β univ T`. -/
theorem gaffkeK_twoPoint [Fintype ι] (γ β : ι → ℝ) (T : Finset ι) :
    gaffkeK (fun i => if i ∈ T then 1 + β i else 1 - γ i) = cubeK γ β Finset.univ T := by
  unfold gaffkeK cubeK
  congr 2
  refine Multiset.map_congr rfl fun i _ => ?_
  simp only [coef]
  split_ifs <;> ring

/-- The posterior's distribution function, `P(L⁺ ≤ x)` for `L⁺ = ∑ U_i ℓ_i + U₀ B`,
`U ~ Dir(1, …, 1)` (Snell and Griffiths), the Dirichlet weights normalized independent unit
exponentials: `P((x − B) E₀ + ∑ (x − ℓ_i) E_i ≥ 0)`. -/
noncomputable def lplusCdf {ι : Type*} [Fintype ι] (B : ℝ) (ℓ : ι → ℝ) (x : ℝ) : ℝ :=
  stat ((x - B) ::ₘ Finset.univ.val.map fun i => x - ℓ i)

/-- **Gaffke's bound** `b*_β`, the least `x` with `P(L⁺ ≤ x) ≥ β` (`LossBound.Gaffke`). -/
noncomputable def gaffkeBound {ι : Type*} [Fintype ι] (B β : ℝ) (ℓ : ι → ℝ) : ℝ :=
  sInf {x | β ≤ lplusCdf B ℓ x}

theorem lplusCdf_mono {ι : Type*} [Fintype ι] (B : ℝ) (ℓ : ι → ℝ) {x y : ℝ} (h : x ≤ y) :
    lplusCdf B ℓ x ≤ lplusCdf B ℓ y := by
  unfold lplusCdf
  calc stat ((x - B) ::ₘ Finset.univ.val.map fun i => x - ℓ i)
      ≤ stat ((y - B) ::ₘ Finset.univ.val.map fun i => x - ℓ i) := stat_cons_le (by linarith) _
    _ ≤ stat ((y - B) ::ₘ Finset.univ.val.map fun i => y - ℓ i) := by
        rw [← Multiset.singleton_add, ← Multiset.singleton_add]
        exact stat_map_mono _ _ fun i _ => by linarith

theorem lplusCdf_nonempty {ι : Type*} [Fintype ι] {B β : ℝ} {ℓ : ι → ℝ} (hℓ : ∀ i, ℓ i ≤ B)
    (hβ1 : β ≤ 1) : ({x | β ≤ lplusCdf B ℓ x} : Set ℝ).Nonempty := by
  refine ⟨B, ?_⟩
  change β ≤ lplusCdf B ℓ B
  rw [lplusCdf, stat_eq_one_of_nonneg]
  · exact hβ1
  · intro r hr
    rcases Multiset.mem_cons.1 hr with rfl | hr
    · simp
    · obtain ⟨i, -, rfl⟩ := Multiset.mem_map.1 hr
      linarith [hℓ i]

/-- **A miss puts `K` low.** If Gaffke's bound falls below `J < B`, then with
`X_i = (B − ℓ_i)/(B − J)`, `K(X) = 1 − P(L⁺ ≤ J) ≤ 1 − β`. -/
theorem gaffkeK_le_of_miss {ι : Type*} [Fintype ι] {B J β : ℝ} {ℓ : ι → ℝ}
    (hℓ : ∀ i, ℓ i ≤ B) (hJB : J < B) (hβ1 : β ≤ 1) (hmiss : gaffkeBound B β ℓ < J) :
    gaffkeK (fun i => (B - ℓ i) / (B - J)) ≤ 1 - β := by
  obtain ⟨x, hx, hxJ⟩ := exists_lt_of_csInf_lt (lplusCdf_nonempty hℓ hβ1) hmiss
  have hcdf : β ≤ lplusCdf B ℓ J := le_trans hx (lplusCdf_mono B ℓ hxJ.le)
  have hBJ : 0 < B - J := sub_pos.2 hJB
  set m : Multiset ℝ := (B - J) ::ₘ Finset.univ.val.map fun i => ℓ i - J with hmdef
  have hK : gaffkeK (fun i => (B - ℓ i) / (B - J)) = stat m := by
    unfold gaffkeK
    rw [← stat_map_mul hBJ, Multiset.map_cons, Multiset.map_map, hmdef]
    congr 2
    · ring
    · exact Multiset.map_congr rfl fun i _ => by
        simp only [Function.comp]
        field_simp
        ring
  have hL : lplusCdf B ℓ J = stat (m.map Neg.neg) := by
    unfold lplusCdf
    rw [hmdef, Multiset.map_cons, Multiset.map_map]
    congr 2
    · ring
    · exact Multiset.map_congr rfl fun i _ => by simp only [Function.comp]; ring
  have hsum := stat_add_stat_neg (m := m) (Multiset.mem_cons_self _ _) hBJ.ne'
  linarith

/-- With every loss at `B`, `L⁺ = B`, and Gaffke's bound is `B`. -/
theorem le_gaffkeBound_of_eq {ι : Type*} [Fintype ι] {B β : ℝ} {ℓ : ι → ℝ} (hℓ : ∀ i, ℓ i = B)
    (hβ0 : 0 < β) (hβ1 : β ≤ 1) : B ≤ gaffkeBound B β ℓ := by
  refine le_csInf (lplusCdf_nonempty (fun i => (hℓ i).le) hβ1) fun x hx => ?_
  by_contra hxB
  have h0 : lplusCdf B ℓ x = 0 := by
    rw [lplusCdf]
    refine stat_eq_zero_of_nonpos (by linarith [not_le.1 hxB]) fun r hr => ?_
    obtain ⟨i, -, rfl⟩ := Multiset.mem_map.1 hr
    linarith [hℓ i, not_le.1 hxB]
  have : β ≤ lplusCdf B ℓ x := hx
  linarith

/-- **The bridge to the cube.** Independent units, each on two values `1 − γ_i` and `1 + β_i`
of a finite index set, are the two-point system `bern` on the high sets. -/
theorem sum_prod_twoPoint {κ : Type*} [Fintype ι] [Fintype κ] (q x : ι → κ → ℝ)
    (γ β : ι → ℝ) (hq1 : ∀ i, ∑ k, q i k = 1)
    (hsupp : ∀ i k, q i k ≠ 0 → x i k = 1 - γ i ∨ x i k = 1 + β i)
    (hhi : ∀ i, ∑ k ∈ Finset.univ.filter (fun k => x i k = 1 + β i), q i k = highChance γ β i)
    (f : (ι → ℝ) → ℝ) :
    ∑ k : ι → κ, (∏ i, q i (k i)) * f (fun i => x i (k i)) =
      ∑ T ∈ (Finset.univ : Finset ι).powerset, bern (highChance γ β) Finset.univ T *
        f (fun i => if i ∈ T then 1 + β i else 1 - γ i) := by
  classical
  set Tk : (ι → κ) → Finset ι := fun k => Finset.univ.filter fun i => x i (k i) = 1 + β i
    with hTk
  rw [← Finset.sum_fiberwise_of_maps_to (s := Finset.univ) (t := Finset.univ.powerset) (g := Tk)
    (fun k _ => Finset.mem_powerset.2 (Finset.subset_univ _))]
  refine Finset.sum_congr rfl fun T _ => ?_
  have hmemT : ∀ k i, (i ∈ Tk k ↔ x i (k i) = 1 + β i) := fun k i => by simp [hTk]
  have hpt : ∀ k ∈ Finset.univ.filter (fun k => Tk k = T),
      (∏ i, q i (k i)) * f (fun i => x i (k i)) =
        (∏ i, q i (k i)) * f (fun i => if i ∈ T then 1 + β i else 1 - γ i) := by
    intro k hk
    have hkT : Tk k = T := (Finset.mem_filter.1 hk).2
    by_cases h0 : ∏ i, q i (k i) = 0
    · rw [h0, zero_mul, zero_mul]
    · congr 2
      funext i
      have hqi : q i (k i) ≠ 0 := fun h => h0 (Finset.prod_eq_zero (Finset.mem_univ i) h)
      have hmem : i ∈ T ↔ x i (k i) = 1 + β i := by rw [← hkT]; exact hmemT k i
      by_cases hiT : i ∈ T
      · rw [if_pos hiT]; exact hmem.1 hiT
      · rw [if_neg hiT]
        rcases hsupp i (k i) hqi with h | h
        · exact h
        · exact absurd (hmem.2 h) hiT
  rw [Finset.sum_congr rfl hpt, ← Finset.sum_mul]
  congr 1
  rw [Finset.sum_filter]
  have hind : ∀ k : ι → κ, (if Tk k = T then ∏ i, q i (k i) else 0) =
      ∏ i, (q i (k i) * if (x i (k i) = 1 + β i ↔ i ∈ T) then 1 else 0) := by
    intro k
    rw [Finset.prod_mul_distrib]
    by_cases hkT : Tk k = T
    · rw [if_pos hkT, Finset.prod_eq_one fun i _ => if_pos (by rw [← hkT]; exact (hmemT k i).symm),
        mul_one]
    · rw [if_neg hkT]
      obtain ⟨i, hi⟩ : ∃ i, ¬ (x i (k i) = 1 + β i ↔ i ∈ T) := by
        by_contra hall
        push Not at hall
        exact hkT (Finset.ext fun i => (hmemT k i).trans (hall i))
      rw [Finset.prod_eq_zero (s := Finset.univ)
        (f := fun i => if (x i (k i) = 1 + β i ↔ i ∈ T) then (1 : ℝ) else 0) (Finset.mem_univ i)
        (if_neg hi), mul_zero]
  rw [Finset.sum_congr rfl fun k _ => hind k,
    ← Fintype.prod_sum (fun i j => q i j * if (x i j = 1 + β i ↔ i ∈ T) then (1 : ℝ) else 0)]
  have hb : bern (highChance γ β) Finset.univ T =
      ∏ i, (if i ∈ T then highChance γ β i else 1 - highChance γ β i) := by
    rw [Finset.prod_ite, bern]
    congr 1
    · congr 1; ext i; simp
    · congr 1; ext i; simp
  rw [hb]
  refine Finset.prod_congr rfl fun i _ => ?_
  by_cases hiT : i ∈ T
  · simp only [hiT, iff_true, if_true]
    rw [← hhi i, Finset.sum_filter]
    exact Finset.sum_congr rfl fun k _ => by split_ifs <;> simp
  · simp only [hiT, iff_false, if_false]
    have hsplit := Finset.sum_filter_add_sum_filter_not Finset.univ
      (fun k => x i k = 1 + β i) (q i)
    rw [hhi i, hq1 i] at hsplit
    rw [show ∑ k, q i k * (if ¬ x i k = 1 + β i then (1 : ℝ) else 0) =
        ∑ k ∈ Finset.univ.filter (fun k => ¬ x i k = 1 + β i), q i k by
      rw [Finset.sum_filter]
      exact Finset.sum_congr rfl fun k _ => by split_ifs <;> simp]
    linarith

omit [DecidableEq ι] in
/-- **The two-point mixture** (Vlassis and Thomas's Lemma 7, for a finite law): a mean-one law on
nonnegative values is a mixture of mean-one laws on two values, one below `1` and one above, and
of the point `1`. The pair `x_a < 1 < x_b` has weight `(x_b − x_a) p_a p_b / M`,
`M = ∑ p (1 − x)⁺`. -/
theorem twoPoint_mixture {κ : Type*} [Fintype κ] (p x : κ → ℝ)
    (hp : ∀ k, 0 ≤ p k) (hx : ∀ k, 0 ≤ x k) (hm : ∑ k, p k * x k = ∑ k, p k) :
    ∃ (w : Option (κ × κ) → ℝ) (q : Option (κ × κ) → κ → ℝ) (γ β : Option (κ × κ) → ℝ),
      (∀ ξ, 0 ≤ w ξ) ∧ (∀ k, p k = ∑ ξ, w ξ * q ξ k) ∧
      ∀ ξ, w ξ ≠ 0 → (∑ k, q ξ k = 1) ∧ 0 ≤ γ ξ ∧ γ ξ ≤ 1 ∧ 0 < β ξ ∧
        (∀ k, q ξ k ≠ 0 → x k = 1 - γ ξ ∨ x k = 1 + β ξ) ∧
        ∑ k ∈ Finset.univ.filter (fun k => x k = 1 + β ξ), q ξ k = γ ξ / (γ ξ + β ξ) := by
  classical
  set M := ∑ k, p k * max 0 (1 - x k) with hM
  have hM' : ∑ k, p k * max 0 (x k - 1) = M := by
    have h : ∀ k, p k * max 0 (x k - 1) - p k * max 0 (1 - x k) = p k * x k - p k := fun k => by
      rcases le_total (x k) 1 with h | h
      · rw [max_eq_left (by linarith), max_eq_right (by linarith)]; ring
      · rw [max_eq_right (by linarith), max_eq_left (by linarith)]; ring
    have := Finset.sum_congr rfl fun k (_ : k ∈ Finset.univ) => h k
    rw [Finset.sum_sub_distrib, Finset.sum_sub_distrib, hm, sub_self] at this
    linarith
  set W := ∑ k, if x k = 1 then p k else 0 with hW
  let w : Option (κ × κ) → ℝ := fun ξ => match ξ with
    | none => W
    | some (a, b) => if x a < 1 ∧ 1 < x b then (x b - x a) * p a * p b / M else 0
  let q : Option (κ × κ) → κ → ℝ := fun ξ k => match ξ with
    | none => if x k = 1 then p k / W else 0
    | some (a, b) => (if k = a then (x b - 1) / (x b - x a) else 0) +
        (if k = b then (1 - x a) / (x b - x a) else 0)
  let γ : Option (κ × κ) → ℝ := fun ξ => match ξ with
    | none => 0
    | some (a, _) => 1 - x a
  let β : Option (κ × κ) → ℝ := fun ξ => match ξ with
    | none => 1
    | some (_, b) => x b - 1
  have hMnn : 0 ≤ M := Finset.sum_nonneg fun k _ => mul_nonneg (hp k) (le_max_left _ _)
  refine ⟨w, q, γ, β, fun ξ => ?_, fun k => ?_, fun ξ hξ => ?_⟩
  · rcases ξ with _ | ⟨a, b⟩
    · exact Finset.sum_nonneg fun k _ => by split_ifs <;> simp [hp k]
    · change 0 ≤ if x a < 1 ∧ 1 < x b then (x b - x a) * p a * p b / M else 0
      split_ifs with h
      · exact div_nonneg (mul_nonneg (mul_nonneg (by linarith [h.1, h.2]) (hp a)) (hp b)) hMnn
      · exact le_rfl
  · -- the decomposition
    rw [Fintype.sum_option, Fintype.sum_prod_type]
    change p k = W * (if x k = 1 then p k / W else 0) + ∑ a, ∑ b,
      (if x a < 1 ∧ 1 < x b then (x b - x a) * p a * p b / M else 0) *
        ((if k = a then (x b - 1) / (x b - x a) else 0) +
          (if k = b then (1 - x a) / (x b - x a) else 0))
    have hsplit : ∀ (u : κ → κ → ℝ) (c₁ c₂ : κ → κ → ℝ),
        ∑ a, ∑ b, u a b * ((if k = a then c₁ a b else 0) + (if k = b then c₂ a b else 0)) =
          ∑ b, u k b * c₁ k b + ∑ a, u a k * c₂ a k := fun u c₁ c₂ => by
      simp only [mul_add, Finset.sum_add_distrib, mul_ite, mul_zero]
      congr 1
      · rw [Finset.sum_comm]
        exact Finset.sum_congr rfl fun b _ => by rw [Finset.sum_ite_eq]; simp
      · exact Finset.sum_congr rfl fun a _ => by rw [Finset.sum_ite_eq]; simp
    rw [hsplit (fun a b => if x a < 1 ∧ 1 < x b then (x b - x a) * p a * p b / M else 0)
      (fun a b => (x b - 1) / (x b - x a)) (fun a b => (1 - x a) / (x b - x a))]
    rcases lt_trichotomy (x k) 1 with hk | hk | hk
    · -- below one: the pairs `(k, b)`
      rw [if_neg hk.ne, mul_zero, zero_add]
      have h2 : ∀ a, (if x a < 1 ∧ 1 < x k then (x k - x a) * p a * p k / M else 0) *
          ((1 - x a) / (x k - x a)) = 0 := fun a => by
        rw [if_neg (by intro h; linarith [h.2])]; simp
      simp only [h2, Finset.sum_const_zero, add_zero]
      have h1 : ∀ b, (if x k < 1 ∧ 1 < x b then (x b - x k) * p k * p b / M else 0) *
          ((x b - 1) / (x b - x k)) = p k / M * (p b * max 0 (x b - 1)) := fun b => by
        by_cases hb : 1 < x b
        · rw [if_pos ⟨hk, hb⟩, max_eq_right (by linarith)]
          field_simp [show x b - x k ≠ 0 by linarith]
        · rw [if_neg (fun h => hb h.2), max_eq_left (by linarith [not_lt.1 hb])]; ring
      simp only [h1, ← Finset.mul_sum, hM']
      by_cases hM0 : M = 0
      · have : p k * max 0 (1 - x k) ≤ M := by
          rw [hM]
          exact Finset.single_le_sum (f := fun j => p j * max 0 (1 - x j))
            (fun j _ => mul_nonneg (hp j) (le_max_left _ _)) (Finset.mem_univ k)
        rw [hM0, max_eq_right (by linarith)] at this
        have hpk : p k = 0 := le_antisymm (by nlinarith [hp k]) (hp k)
        rw [hpk, zero_div, zero_mul]
      · field_simp
    · -- at one: the point `1`
      rw [if_pos hk]
      have h1 : ∀ b, (if x k < 1 ∧ 1 < x b then (x b - x k) * p k * p b / M else 0) *
          ((x b - 1) / (x b - x k)) = 0 := fun b => by
        rw [if_neg (by intro h; linarith [h.1])]; simp
      have h2 : ∀ a, (if x a < 1 ∧ 1 < x k then (x k - x a) * p a * p k / M else 0) *
          ((1 - x a) / (x k - x a)) = 0 := fun a => by
        rw [if_neg (by intro h; linarith [h.2])]; simp
      simp only [h1, h2, Finset.sum_const_zero, add_zero]
      by_cases hW0 : W = 0
      · have : p k ≤ W := by
          rw [hW]
          calc p k = if x k = 1 then p k else 0 := by rw [if_pos hk]
            _ ≤ ∑ j, if x j = 1 then p j else 0 :=
              Finset.single_le_sum (f := fun j => if x j = 1 then p j else 0)
                (fun j _ => by split_ifs <;> simp [hp j]) (Finset.mem_univ k)
        rw [hW0] at this
        rw [hW0, zero_mul]
        exact le_antisymm this (hp k)
      · field_simp
    · -- above one: the pairs `(a, k)`
      rw [if_neg hk.ne', mul_zero, zero_add]
      have h1 : ∀ b, (if x k < 1 ∧ 1 < x b then (x b - x k) * p k * p b / M else 0) *
          ((x b - 1) / (x b - x k)) = 0 := fun b => by
        rw [if_neg (by intro h; linarith [h.1])]; simp
      simp only [h1, Finset.sum_const_zero, zero_add]
      have h2 : ∀ a, (if x a < 1 ∧ 1 < x k then (x k - x a) * p a * p k / M else 0) *
          ((1 - x a) / (x k - x a)) = p k / M * (p a * max 0 (1 - x a)) := fun a => by
        by_cases ha : x a < 1
        · rw [if_pos ⟨ha, hk⟩, max_eq_right (by linarith)]
          field_simp [show x k - x a ≠ 0 by linarith]
        · rw [if_neg (fun h => ha h.1), max_eq_left (by linarith [not_lt.1 ha])]; ring
      simp only [h2, ← Finset.mul_sum, ← hM]
      by_cases hM0 : M = 0
      · have : p k * max 0 (x k - 1) ≤ M := by
          rw [← hM']
          exact Finset.single_le_sum (f := fun j => p j * max 0 (x j - 1))
            (fun j _ => mul_nonneg (hp j) (le_max_left _ _)) (Finset.mem_univ k)
        rw [hM0, max_eq_right (by linarith)] at this
        have hpk : p k = 0 := le_antisymm (by nlinarith [hp k]) (hp k)
        rw [hpk, zero_div, zero_mul]
      · field_simp
  · -- each component with weight is a mean-one two-point law
    rcases ξ with _ | ⟨a, b⟩
    · have hW0 : W ≠ 0 := hξ
      refine ⟨?_, le_rfl, zero_le_one, one_pos, fun k hk => ?_, ?_⟩
      · change ∑ k, (if x k = 1 then p k / W else 0) = 1
        rw [show (fun k => if x k = 1 then p k / W else 0) = fun k =>
          (if x k = 1 then p k else 0) / W from funext fun k => by split_ifs <;> simp,
          ← Finset.sum_div, ← hW, div_self hW0]
      · left
        change x k = 1 - 0
        by_contra h
        exact hk (show (if x k = 1 then p k / W else 0) = 0 from if_neg (by simpa using h))
      · change ∑ k ∈ Finset.univ.filter (fun k => x k = 1 + 1), (if x k = 1 then p k / W else 0) =
          0 / (0 + 1)
        rw [zero_div]
        exact Finset.sum_eq_zero fun k hk => if_neg (by
          have := (Finset.mem_filter.1 hk).2; linarith)
    · have hab : x a < 1 ∧ 1 < x b := by
        by_contra h
        exact hξ (show (if x a < 1 ∧ 1 < x b then (x b - x a) * p a * p b / M else 0) = 0 from
          if_neg h)
      have hD : x b - x a ≠ 0 := by linarith [hab.1, hab.2]
      have hne : a ≠ b := fun h => by rw [h] at hab; linarith [hab.1, hab.2]
      refine ⟨?_, by change 0 ≤ 1 - x a; linarith [hab.1], by change 1 - x a ≤ 1; linarith [hx a],
        by change 0 < x b - 1; linarith [hab.2], fun k hk => ?_, ?_⟩
      · change ∑ k, ((if k = a then (x b - 1) / (x b - x a) else 0) +
          (if k = b then (1 - x a) / (x b - x a) else 0)) = 1
        rw [Finset.sum_add_distrib, Finset.sum_ite_eq', Finset.sum_ite_eq',
          if_pos (Finset.mem_univ _), if_pos (Finset.mem_univ _)]
        field_simp
        ring
      · have hk' : (if k = a then (x b - 1) / (x b - x a) else 0) +
            (if k = b then (1 - x a) / (x b - x a) else 0) ≠ 0 := hk
        by_cases hka : k = a
        · left; change x k = 1 - (1 - x a); rw [hka]; ring
        · right; change x k = 1 + (x b - 1)
          by_cases hkb : k = b
          · rw [hkb]; ring
          · exact absurd (by rw [if_neg hka, if_neg hkb, add_zero]) hk'
      · change ∑ k ∈ Finset.univ.filter (fun k => x k = 1 + (x b - 1)),
          ((if k = a then (x b - 1) / (x b - x a) else 0) +
            (if k = b then (1 - x a) / (x b - x a) else 0)) = (1 - x a) / (1 - x a + (x b - 1))
        rw [Finset.sum_add_distrib, Finset.sum_ite_eq', Finset.sum_ite_eq',
          if_neg (by simp only [Finset.mem_filter, Finset.mem_univ, true_and]; intro h; linarith),
          if_pos (by simp only [Finset.mem_filter, Finset.mem_univ, true_and]; ring), zero_add]
        congr 1
        ring

/-- **Vlassis and Thomas's Theorem 1, for finite laws.** Independent units with nonnegative
values and means at most `1`: `K(X) ≤ α` with chance at most `α` (`max α 0` for any real `α`).
Rescale to mean one (`K` falls as the values rise), split each law into two-point laws
(`twoPoint_mixture`), and each product of two-point laws is `twoPoint_valid`. -/
theorem gaffkeK_valid {κ : Type*} [Fintype ι] [Fintype κ]
    (p x : ι → κ → ℝ) (hp : ∀ i k, 0 ≤ p i k) (hp1 : ∀ i, ∑ k, p i k = 1)
    (hx : ∀ i k, 0 ≤ x i k) (hmean : ∀ i, ∑ k, p i k * x i k ≤ 1) (α : ℝ) :
    ∑ k : ι → κ, (∏ i, p i (k i)) * (if gaffkeK (fun i => x i (k i)) ≤ α then 1 else 0) ≤
      max α 0 := by
  classical
  -- rescale to mean one
  set mu : ι → ℝ := fun i => ∑ k, p i k * x i k with hmu
  set x' : ι → κ → ℝ := fun i k => if mu i = 0 then 1 + x i k else x i k / mu i with hx'
  have hmu0 : ∀ i, 0 ≤ mu i := fun i => Finset.sum_nonneg fun k _ => mul_nonneg (hp i k) (hx i k)
  have hxx' : ∀ i k, x i k ≤ x' i k := fun i k => by
    simp only [hx']
    split_ifs with h
    · linarith
    · rw [le_div_iff₀ (lt_of_le_of_ne (hmu0 i) (Ne.symm h))]
      nlinarith [hx i k, hmean i]
  have hx'0 : ∀ i k, 0 ≤ x' i k := fun i k => (hx i k).trans (hxx' i k)
  have hm' : ∀ i, ∑ k, p i k * x' i k = ∑ k, p i k := fun i => by
    rw [hp1 i]
    simp only [hx']
    split_ifs with h
    · simp only [mul_add, mul_one, Finset.sum_add_distrib, hp1 i]
      have : mu i = 0 := h
      rw [hmu] at this
      linarith
    · simp only [mul_div_assoc', ← Finset.sum_div]
      exact div_self h
  calc ∑ k : ι → κ, (∏ i, p i (k i)) * (if gaffkeK (fun i => x i (k i)) ≤ α then 1 else 0)
      ≤ ∑ k : ι → κ, (∏ i, p i (k i)) * (if gaffkeK (fun i => x' i (k i)) ≤ α then 1 else 0) :=
        Finset.sum_le_sum fun k _ => mul_le_mul_of_nonneg_left (by
          split_ifs with h1 h2
          · exact le_rfl
          · exact absurd ((gaffkeK_anti fun i => hxx' i (k i)).trans h1) h2
          · norm_num
          · exact le_rfl) (Finset.prod_nonneg fun i _ => hp i (k i))
    _ ≤ max α 0 := by
      -- split every unit's law into two-point laws
      choose w q γ β hw hdec hcomp using fun i => twoPoint_mixture (p i) (x' i) (hp i) (hx'0 i)
        (hm' i)
      have hwsum : ∀ i, ∑ ξ, w i ξ = 1 := fun i => by
        rw [← hp1 i]
        calc ∑ ξ, w i ξ = ∑ ξ, w i ξ * ∑ k, q i ξ k := Finset.sum_congr rfl fun ξ _ => by
              by_cases h : w i ξ = 0
              · rw [h, zero_mul]
              · rw [(hcomp i ξ h).1, mul_one]
          _ = ∑ k, p i k := by
              simp only [Finset.mul_sum]
              rw [Finset.sum_comm]
              exact Finset.sum_congr rfl fun k _ => (hdec i k).symm
      set F : (ι → κ) → ℝ := fun k => if gaffkeK (fun i => x' i (k i)) ≤ α then 1 else 0 with hF
      have hmix : ∑ k : ι → κ, (∏ i, p i (k i)) * F k =
          ∑ ξ : ι → Option (κ × κ), (∏ i, w i (ξ i)) *
            ∑ k : ι → κ, (∏ i, q i (ξ i) (k i)) * F k := by
        calc ∑ k : ι → κ, (∏ i, p i (k i)) * F k
            = ∑ k : ι → κ, (∑ ξ : ι → Option (κ × κ), ∏ i, (w i (ξ i) * q i (ξ i) (k i))) * F k :=
              Finset.sum_congr rfl fun k _ => by
                rw [← Fintype.prod_sum (fun i ξ => w i ξ * q i ξ (k i))]
                exact congrArg (· * F k) (Finset.prod_congr rfl fun i _ => hdec i (k i))
          _ = _ := by
              simp only [Finset.sum_mul, Finset.mul_sum, Finset.prod_mul_distrib]
              rw [Finset.sum_comm]
              exact Finset.sum_congr rfl fun ξ _ => Finset.sum_congr rfl fun k _ => by ring
      rw [hmix]
      calc ∑ ξ : ι → Option (κ × κ), (∏ i, w i (ξ i)) *
            ∑ k : ι → κ, (∏ i, q i (ξ i) (k i)) * F k
          ≤ ∑ ξ : ι → Option (κ × κ), (∏ i, w i (ξ i)) * max α 0 := by
            refine Finset.sum_le_sum fun ξ _ => ?_
            by_cases h0 : ∏ i, w i (ξ i) = 0
            · rw [h0, zero_mul, zero_mul]
            · have hwi : ∀ i, w i (ξ i) ≠ 0 := fun i h =>
                h0 (Finset.prod_eq_zero (Finset.mem_univ i) h)
              refine mul_le_mul_of_nonneg_left ?_ (Finset.prod_nonneg fun i _ => hw i (ξ i))
              have hbr := sum_prod_twoPoint (fun i => q i (ξ i)) x' (fun i => γ i (ξ i))
                (fun i => β i (ξ i)) (fun i => (hcomp i (ξ i) (hwi i)).1)
                (fun i => (hcomp i (ξ i) (hwi i)).2.2.2.2.1)
                (fun i => (hcomp i (ξ i) (hwi i)).2.2.2.2.2)
                (fun y => if gaffkeK y ≤ α then (1 : ℝ) else 0)
              refine le_of_eq_of_le hbr ?_
              simp only [gaffkeK_twoPoint]
              exact twoPoint_valid (fun i => (hcomp i (ξ i) (hwi i)).2.1)
                (fun i => (hcomp i (ξ i) (hwi i)).2.2.1)
                (fun i => (hcomp i (ξ i) (hwi i)).2.2.2.1) α
        _ = max α 0 := by
            rw [← Finset.sum_mul, ← Fintype.prod_sum]
            simp [hwsum]

/-- **Gaffke's bound misses with chance at most `1 − β`** (Vlassis and Thomas's Theorem 1, through
`P(L⁺ ≤ J) = 1 − K(X)` with `X_i = (B − ℓ_i)/(B − J)`): independent units, each with its own loss
law on finitely many values at most `B` and mean at least `J ≤ B`; `J ≤ B` matters only with no
units. -/
theorem gaffke_miss_le {κ : Type*} [Fintype ι] [Fintype κ]
    (p v : ι → κ → ℝ) (hp : ∀ i k, 0 ≤ p i k) (hp1 : ∀ i, ∑ k, p i k = 1) {B J β : ℝ}
    (hv : ∀ i k, v i k ≤ B) (hJ : ∀ i, J ≤ ∑ k, p i k * v i k) (hJB0 : J ≤ B) (hβ1 : β ≤ 1) :
    ∑ k : ι → κ, (∏ i, p i (k i)) *
      (if gaffkeBound B β (fun i => v i (k i)) < J then 1 else 0) ≤ 1 - β := by
  classical
  have hprod : ∀ k : ι → κ, 0 ≤ ∏ i, p i (k i) := fun k => Finset.prod_nonneg fun i _ => hp i (k i)
  have htot : ∑ k : ι → κ, ∏ i, p i (k i) = 1 := by
    rw [← Fintype.prod_sum]; simp [hp1]
  rcases le_or_gt β 0 with hβ0 | hβ0
  · calc ∑ k : ι → κ, (∏ i, p i (k i)) *
          (if gaffkeBound B β (fun i => v i (k i)) < J then 1 else 0)
        ≤ ∑ k : ι → κ, ∏ i, p i (k i) := Finset.sum_le_sum fun k _ => by
          split_ifs <;> simp [hprod k]
      _ ≤ 1 - β := by rw [htot]; linarith
  have hmeanB : ∀ i, ∑ k, p i k * v i k ≤ B := fun i => by
    calc ∑ k, p i k * v i k ≤ ∑ k, p i k * B :=
          Finset.sum_le_sum fun k _ => mul_le_mul_of_nonneg_left (hv i k) (hp i k)
      _ = B := by rw [← Finset.sum_mul, hp1 i, one_mul]
  rcases lt_or_ge J B with hJB | hJB
  · -- `X_i = (B − ℓ_i)/(B − J)`, mean at most one, and Theorem 1 at `α = 1 − β`
    have hBJ : 0 < B - J := sub_pos.2 hJB
    have hval := gaffkeK_valid p (fun i k => (B - v i k) / (B - J)) hp hp1
      (fun i k => div_nonneg (sub_nonneg.2 (hv i k)) hBJ.le) (fun i => by
        simp only [mul_div_assoc', ← Finset.sum_div]
        rw [div_le_one hBJ]
        simp only [mul_sub, Finset.sum_sub_distrib, ← Finset.sum_mul, hp1 i, one_mul]
        linarith [hJ i]) (1 - β)
    rw [max_eq_left (by linarith)] at hval
    refine le_trans (Finset.sum_le_sum fun k _ => mul_le_mul_of_nonneg_left ?_ (hprod k)) hval
    split_ifs with h1 h2
    · exact le_rfl
    · exact absurd (gaffkeK_le_of_miss (fun i => hv i (k i)) hJB hβ1 h1) h2
    · norm_num
    · exact le_rfl
  · -- `J = B`: every loss with chance is `B`, and the bound is `B`
    refine le_trans (le_of_eq (Finset.sum_eq_zero fun k _ => ?_)) (by linarith)
    by_cases h0 : ∏ i, p i (k i) = 0
    · rw [h0, zero_mul]
    · have hall : ∀ i, v i (k i) = B := fun i => by
        have hpi : p i (k i) ≠ 0 := fun h => h0 (Finset.prod_eq_zero (Finset.mem_univ i) h)
        have hsum : ∑ k', p i k' * (B - v i k') = 0 := by
          have := hJ i
          have := hmeanB i
          simp only [mul_sub, Finset.sum_sub_distrib, ← Finset.sum_mul, hp1 i, one_mul]
          linarith
        have := (Finset.sum_eq_zero_iff_of_nonneg fun k' _ =>
          mul_nonneg (hp i k') (sub_nonneg.2 (hv i k'))).1 hsum (k i) (Finset.mem_univ _)
        rcases mul_eq_zero.1 this with h | h
        · exact absurd h hpi
        · linarith
      rw [if_neg (not_lt.2 ((le_of_eq (le_antisymm hJB0 hJB)).trans
        (le_gaffkeBound_of_eq hall hβ0 hβ1))), mul_zero]

/-- **The `β`-bound's calibration at `n` units**: `P(J ≤ b*_β) ≥ β` (Gaffke's conjecture, Vlassis
and Thomas's Theorem 1), for independent units with their own finite loss laws at most `B`, each
mean at least `J ≤ B`. `gaffke_one` is the case of one unit. -/
theorem gaffke_calibration {κ : Type*} [Fintype ι] [Fintype κ]
    (p v : ι → κ → ℝ) (hp : ∀ i k, 0 ≤ p i k) (hp1 : ∀ i, ∑ k, p i k = 1) {B J β : ℝ}
    (hv : ∀ i k, v i k ≤ B) (hJ : ∀ i, J ≤ ∑ k, p i k * v i k) (hJB0 : J ≤ B) (hβ1 : β ≤ 1) :
    β ≤ ∑ k : ι → κ, (∏ i, p i (k i)) *
      (if J ≤ gaffkeBound B β (fun i => v i (k i)) then 1 else 0) := by
  classical
  have htot : ∑ k : ι → κ, ∏ i, p i (k i) = 1 := by
    rw [← Fintype.prod_sum]; simp [hp1]
  have hpt : ∀ k : ι → κ, (∏ i, p i (k i)) *
      (if J ≤ gaffkeBound B β (fun i => v i (k i)) then 1 else 0) + (∏ i, p i (k i)) *
        (if gaffkeBound B β (fun i => v i (k i)) < J then 1 else 0) = ∏ i, p i (k i) := fun k => by
    split_ifs with h1 h2 <;> first | (exfalso; linarith) | ring
  have hcompl := Finset.sum_congr rfl fun k (_ : k ∈ (Finset.univ : Finset (ι → κ))) => hpt k
  rw [Finset.sum_add_distrib, htot] at hcompl
  have := gaffke_miss_le p v hp hp1 hv hJ hJB0 hβ1
  linarith

/-- Gaffke's bound rises with the losses. -/
theorem gaffkeBound_mono {ι : Type*} [Fintype ι] {B β : ℝ} (hβ0 : 0 < β) (hβ1 : β ≤ 1)
    {ℓ ℓ' : ι → ℝ} (hℓ' : ∀ i, ℓ' i ≤ B) (h : ∀ i, ℓ i ≤ ℓ' i) :
    gaffkeBound B β ℓ ≤ gaffkeBound B β ℓ' := by
  have hℓ : ∀ i, ℓ i ≤ B := fun i => (h i).trans (hℓ' i)
  refine csInf_le_csInf ⟨B - ∑ i, (B - ℓ i) - 1, fun x hx => ?_⟩
    (lplusCdf_nonempty hℓ' hβ1) fun x hx => ?_
  · by_contra hlt
    rw [not_le] at hlt
    have hsingle : ∀ i, B - ℓ i ≤ ∑ j, (B - ℓ j) := fun i =>
      Finset.single_le_sum (f := fun j => B - ℓ j) (fun j _ => sub_nonneg.2 (hℓ j))
        (Finset.mem_univ i)
    have hsum0 : 0 ≤ ∑ j, (B - ℓ j) := Finset.sum_nonneg fun j _ => sub_nonneg.2 (hℓ j)
    have h0 : lplusCdf B ℓ x = 0 := by
      rw [lplusCdf]
      refine stat_eq_zero_of_nonpos (by linarith) fun r hr => ?_
      obtain ⟨i, -, rfl⟩ := Multiset.mem_map.1 hr
      linarith [hsingle i]
    have : β ≤ lplusCdf B ℓ x := hx
    linarith
  · change β ≤ lplusCdf B ℓ x
    refine le_trans hx ?_
    unfold lplusCdf
    rw [← Multiset.singleton_add, ← Multiset.singleton_add]
    exact stat_map_mono _ _ fun i _ => by linarith [h i]

/-- **The high-probability rule** (Snell and Griffiths's `λ^β_hpd`, `RiskControl.HighProbability`):
the first threshold on the grid whose Gaffke bound is at most `α`, the last when none is. -/
noncomputable def hpdRule {K : ℕ} {ι : Type*} [Fintype ι] (B α β : ℝ)
    (ℓ : ι → Fin (K + 1) → ℝ) : Fin (K + 1) :=
  open Classical in firstFeasible fun l => gaffkeBound B β (fun i => ℓ i l) ≤ α

/-- **What `RiskControl.HighProbability` reports.** Independent calibration units, each with its
own loss law on finitely many loss curves, nonincreasing in the threshold and at most `B`; `J l` at
most `B` and at most every unit's mean loss at threshold `l` (the risk, for exchangeable units).
The rule meets `α` at a threshold whose risk exceeds `α` with chance at most `1 − β`: such a miss
puts the fixed threshold `l₀`, the last with risk above `α`, below its risk (`fail_le_fixed`), and
that has chance at most `1 − β` (`gaffke_miss_le`). -/
theorem highProbability_risk {κ : Type*} [Fintype ι] [Fintype κ] {K : ℕ}
    (p : ι → κ → ℝ) (hp : ∀ i k, 0 ≤ p i k) (hp1 : ∀ i, ∑ k, p i k = 1)
    (v : ι → κ → Fin (K + 1) → ℝ) {B α β : ℝ} (hv : ∀ i k l, v i k l ≤ B)
    (hanti : ∀ i k, Antitone (v i k)) (J : Fin (K + 1) → ℝ)
    (hJ : ∀ i l, J l ≤ ∑ k, p i k * v i k l) (hJB : ∀ l, J l ≤ B) (hβ1 : β ≤ 1) :
    ∑ k : ι → κ, (∏ i, p i (k i)) *
      (if gaffkeBound B β (fun i => v i (k i) (hpdRule B α β fun i => v i (k i))) ≤ α ∧
          α < J (hpdRule B α β fun i => v i (k i)) then 1 else 0) ≤ 1 - β := by
  classical
  have hprod : ∀ k : ι → κ, 0 ≤ ∏ i, p i (k i) := fun k => Finset.prod_nonneg fun i _ => hp i (k i)
  have htot : ∑ k : ι → κ, ∏ i, p i (k i) = 1 := by
    rw [← Fintype.prod_sum]; simp [hp1]
  rcases le_or_gt β 0 with hβ0 | hβ0
  · calc _ ≤ ∑ k : ι → κ, ∏ i, p i (k i) := Finset.sum_le_sum fun k _ => by
          split_ifs <;> simp [hprod k]
      _ ≤ 1 - β := by rw [htot]; linarith
  by_cases hex : ∃ l, α < J l
  · set S := Finset.univ.filter fun l => α < J l with hS
    have hSne : S.Nonempty := let ⟨l, hl⟩ := hex; ⟨l, Finset.mem_filter.2 ⟨Finset.mem_univ _, hl⟩⟩
    set l₀ := S.max' hSne with hl₀def
    have hJ₀ : α < J l₀ := (Finset.mem_filter.1 (S.max'_mem hSne)).2
    have hl₀ : ∀ l', α < J l' → l' ≤ l₀ := fun l' h =>
      S.le_max' l' (Finset.mem_filter.2 ⟨Finset.mem_univ _, h⟩)
    refine le_trans (Finset.sum_le_sum fun k _ => mul_le_mul_of_nonneg_left ?_ (hprod k))
      (gaffke_miss_le p (fun i k => v i k l₀) hp hp1 (fun i k => hv i k l₀) (fun i => hJ i l₀)
        (hJB l₀) hβ1)
    have hb : Antitone fun l => gaffkeBound B β (fun i => v i (k i) l) := fun l l' hll' =>
      gaffkeBound_mono hβ0 hβ1 (fun i => hv i (k i) l) fun i => hanti i (k i) hll'
    split_ifs with h1 h2
    · exact le_rfl
    · exact absurd (fail_le_fixed hb h1.1 hJ₀ hl₀ h1.2) h2
    · norm_num
    · exact le_rfl
  · push Not at hex
    refine le_trans (le_of_eq (Finset.sum_eq_zero fun k _ => ?_)) (by linarith)
    rw [if_neg (fun h => absurd h.2 (not_lt.2 (hex _))), mul_zero]

end VlassisThomas

section ProperScores

variable {κ : Type*} [Fintype κ]

/-- Gibbs' inequality, one outcome at a time: `p log q − p log p ≤ q − p`, for `p, q ≥ 0` with
`q > 0` where `p > 0`. -/
theorem gibbs_term {p q : ℝ} (hp : 0 ≤ p) (hq : 0 ≤ q) (hpq : 0 < p → 0 < q) :
    p * Real.log q - p * Real.log p ≤ q - p := by
  rcases hp.lt_or_eq with h | h
  · have hq' := hpq h
    have hlog := Real.log_le_sub_one_of_pos (div_pos hq' h)
    rw [Real.log_div hq'.ne' h.ne'] at hlog
    have h2 := mul_le_mul_of_nonneg_left hlog h.le
    have h3 : p * (q / p - 1) = q - p := by field_simp
    rw [h3] at h2
    linarith
  · rw [← h]
    simpa using hq

/-- Strictly, where `p ≠ q`. -/
theorem gibbs_term_lt {p q : ℝ} (hp : 0 ≤ p) (hq : 0 ≤ q) (hpq : 0 < p → 0 < q) (hne : p ≠ q) :
    p * Real.log q - p * Real.log p < q - p := by
  rcases hp.lt_or_eq with h | h
  · have hq' := hpq h
    have hne1 : q / p ≠ 1 := fun h1 => hne ((div_eq_one_iff_eq h.ne').1 h1).symm
    have hlog := Real.log_lt_sub_one_of_pos (div_pos hq' h) hne1
    rw [Real.log_div hq'.ne' h.ne'] at hlog
    have h2 := mul_lt_mul_of_pos_left hlog h
    have h3 : p * (q / p - 1) = q - p := by field_simp
    rw [h3] at h2
    linarith
  · subst h
    have : 0 < q := lt_of_le_of_ne hq hne
    simpa using this

/-- **The log-loss is proper** (Gressmann et al.'s Proposition 2.5 (i), through Gibbs' inequality):
for a law `p` on finitely many outcomes and a prediction `q` (nonnegative, total at most one,
positive where `p` is), the expected log-loss of `q` is at least that of `p`. -/
theorem logLoss_proper (p q : κ → ℝ) (hp : ∀ y, 0 ≤ p y) (hp1 : ∑ y, p y = 1)
    (hq : ∀ y, 0 ≤ q y) (hq1 : ∑ y, q y ≤ 1) (hpq : ∀ y, 0 < p y → 0 < q y) :
    ∑ y, p y * -Real.log (p y) ≤ ∑ y, p y * -Real.log (q y) := by
  have h := sum_le_sum fun y (_ : y ∈ univ) => gibbs_term (hp y) (hq y) (hpq y)
  rw [sum_sub_distrib, sum_sub_distrib, hp1] at h
  simp only [mul_neg, sum_neg_distrib]
  linarith

/-- **Strictly proper**: for `q` a law other than `p`, strictly. -/
theorem logLoss_strict (p q : κ → ℝ) (hp : ∀ y, 0 ≤ p y) (hp1 : ∑ y, p y = 1)
    (hq : ∀ y, 0 ≤ q y) (hq1 : ∑ y, q y = 1) (hpq : ∀ y, 0 < p y → 0 < q y) (hne : p ≠ q) :
    ∑ y, p y * -Real.log (p y) < ∑ y, p y * -Real.log (q y) := by
  obtain ⟨y₀, hy₀⟩ := Function.ne_iff.1 hne
  have h := sum_lt_sum (fun y (_ : y ∈ univ) => gibbs_term (hp y) (hq y) (hpq y))
    ⟨y₀, mem_univ _, gibbs_term_lt (hp y₀) (hq y₀) (hpq y₀) hy₀⟩
  rw [sum_sub_distrib, sum_sub_distrib, hp1, hq1] at h
  simp only [mul_neg, sum_neg_distrib]
  linarith

/-- **The integrated squared loss's gap** (Gressmann et al.'s Lemma 2.4 (ii); the loss is minus
Gneiting and Raftery's quadratic score, and the gap is its Bregman divergence): with
`L(q, y) = −2 q(y) + ‖q‖²`, the multi-class Brier loss, the expected loss of `q` exceeds that of
the law `p` by `Σ (p − q)²`. -/
theorem gneitingLoss_gap (p q : κ → ℝ) (hp1 : ∑ y, p y = 1) :
    ∑ y, p y * (-2 * q y + ∑ j, q j ^ 2) - ∑ y, p y * (-2 * p y + ∑ j, p j ^ 2) =
      ∑ y, (p y - q y) ^ 2 := by
  have key : ∀ r : κ → ℝ,
      ∑ y, p y * (-2 * r y + ∑ j, r j ^ 2) = -2 * ∑ y, p y * r y + ∑ j, r j ^ 2 := by
    intro r
    rw [show (∑ y, p y * (-2 * r y + ∑ j, r j ^ 2)) =
        ∑ y, (-2 * (p y * r y) + p y * ∑ j, r j ^ 2) from sum_congr rfl fun y _ => by ring,
      sum_add_distrib, ← mul_sum, ← sum_mul, hp1, one_mul]
  rw [key, key]
  have h1 : ∑ y, (p y - q y) ^ 2 = ∑ y, p y * p y - 2 * ∑ y, p y * q y + ∑ y, q y ^ 2 := by
    rw [mul_sum, ← sum_sub_distrib, ← sum_add_distrib]
    exact sum_congr rfl fun y _ => by ring
  have h2 : ∑ y, p y ^ 2 = ∑ y, p y * p y := sum_congr rfl fun y _ => by ring
  rw [h1, h2]
  ring

/-- **Strictly proper** (Gressmann et al.'s Proposition 2.5 (ii)), on finitely many outcomes. -/
theorem gneitingLoss_strict (p q : κ → ℝ) (hp1 : ∑ y, p y = 1) (hne : p ≠ q) :
    ∑ y, p y * (-2 * p y + ∑ j, p j ^ 2) < ∑ y, p y * (-2 * q y + ∑ j, q j ^ 2) := by
  obtain ⟨y₀, hy₀⟩ := Function.ne_iff.1 hne
  have hpos : 0 < ∑ y, (p y - q y) ^ 2 :=
    sum_pos' (fun y _ => sq_nonneg _)
      ⟨y₀, mem_univ _, lt_of_le_of_ne (sq_nonneg _) (pow_ne_zero 2 (sub_ne_zero.2 hy₀)).symm⟩
  have := gneitingLoss_gap p q hp1
  linarith

/-- **Proper scores rank conditional predictions** (Gressmann et al.'s Proposition 2.7): a loss
proper for every law of a class has the true conditional laws minimizing the expected loss over
features too. -/
theorem proper_conditional {χ : Type*} [Fintype χ] (L : (κ → ℝ) → κ → ℝ)
    (law : (κ → ℝ) → Prop) (hL : ∀ p q, law p → law q → ∑ y, p y * L p y ≤ ∑ y, p y * L q y)
    (w : χ → ℝ) (hw : ∀ x, 0 ≤ w x) (P R : χ → κ → ℝ) (hP : ∀ x, law (P x))
    (hR : ∀ x, law (R x)) :
    ∑ x, w x * ∑ y, P x y * L (P x) y ≤ ∑ x, w x * ∑ y, P x y * L (R x) y :=
  sum_le_sum fun x _ => mul_le_mul_of_nonneg_left (hL _ _ (hP x) (hR x)) (hw x)

end ProperScores

/-! ### RoBAS-Full's conformal sets are intervals for `n ≥ 4`

Commissioned 2026-11-14: RoBAS-Full's sets are intervals from `n = 4`, so that
`ResidualShrinkage.Intervals` can serve it as it serves RoBAS-EB. Sources: RoBAS
(arXiv 2607.04236v2), its appendix's RoBAS-Full derivation and Bersson and Hoff's Lemmas 1–2.
Unit `j` is accepted where `(r_j − z) β_j(z) ≥ 0`. `β_j` rises through each of its zeros, which
reduces to an inequality in Dawson's function; its four bounds (`RobasDawsonBounds`) are certified
by interval arithmetic in `conformal-checks.py`, which reads the constants from this file. Below
`n = 4` the inequality fails, and the sets need not be intervals. -/

section RobasFull

open Real Set Filter MeasureTheory Topology

/-- `K_j(x) = ∫₀¹ v^{2j} e^{x v²} dv`. -/
noncomputable def kummerMoment (j : ℕ) (x : ℝ) : ℝ :=
  ∫ v in (0 : ℝ)..1, v ^ (2 * j) * exp (x * v ^ 2)

theorem kummerMoment_continuous_integrand (j : ℕ) (x : ℝ) :
    Continuous fun v : ℝ => v ^ (2 * j) * exp (x * v ^ 2) := by fun_prop

theorem kummerMoment_pos (j : ℕ) (x : ℝ) : 0 < kummerMoment j x :=
  intervalIntegral.intervalIntegral_pos_of_pos_on
    ((kummerMoment_continuous_integrand j x).intervalIntegrable _ _)
    (fun _ hv => mul_pos (pow_pos hv.1 _) (exp_pos _)) zero_lt_one

theorem kummerMoment_succ_le (j : ℕ) (x : ℝ) : kummerMoment (j + 1) x ≤ kummerMoment j x := by
  refine intervalIntegral.integral_mono_on zero_le_one
    ((kummerMoment_continuous_integrand _ x).intervalIntegrable _ _)
    ((kummerMoment_continuous_integrand j x).intervalIntegrable _ _) fun v hv => ?_
  have h1 : v ^ (2 * (j + 1)) ≤ v ^ (2 * j) := by
    rw [show 2 * (j + 1) = 2 * j + 2 by ring, pow_add]
    exact mul_le_of_le_one_right (pow_nonneg hv.1 _) (pow_le_one₀ hv.1 hv.2)
  exact mul_le_mul_of_nonneg_right h1 (exp_pos _).le

theorem hasDerivAt_kummerMoment (j : ℕ) (x : ℝ) :
    HasDerivAt (kummerMoment j) (kummerMoment (j + 1) x) x := by
  have hc : ∀ y : ℝ, IntegrableOn (fun v : ℝ => v ^ (2 * j) * exp (y * v ^ 2)) (Ioc 0 1) := fun y =>
    ((kummerMoment_continuous_integrand j y).integrableOn_Icc).mono_set Ioc_subset_Icc_self
  have key := hasDerivAt_integral_of_dominated_loc_of_deriv_le
    (μ := volume.restrict (Ioc (0 : ℝ) 1))
    (F := fun y v => v ^ (2 * j) * exp (y * v ^ 2))
    (F' := fun y v => v ^ (2 * (j + 1)) * exp (y * v ^ 2)) (x₀ := x)
    (bound := fun _ => exp (|x| + 1)) (s := Metric.ball x 1) (Metric.ball_mem_nhds x one_pos)
    (Eventually.of_forall fun y => (kummerMoment_continuous_integrand j y).aestronglyMeasurable)
    (hc x) (kummerMoment_continuous_integrand (j + 1) x).aestronglyMeasurable
    ((ae_restrict_iff' measurableSet_Ioc).2 (Eventually.of_forall fun v hv y hy => by
      have hy' : |y| < |x| + 1 := by
        have := Metric.mem_ball.1 hy
        rw [Real.dist_eq] at this
        linarith [abs_sub_abs_le_abs_sub y x]
      have hv2 : v ^ 2 ≤ 1 := pow_le_one₀ hv.1.le hv.2
      rw [Real.norm_eq_abs, abs_mul, abs_of_nonneg (exp_pos _).le,
        abs_of_nonneg (pow_nonneg hv.1.le _)]
      calc v ^ (2 * (j + 1)) * exp (y * v ^ 2) ≤ 1 * exp (|x| + 1) := by
            refine mul_le_mul (pow_le_one₀ hv.1.le hv.2) (exp_le_exp.2 ?_) (exp_pos _).le
              zero_le_one
            calc y * v ^ 2 ≤ |y| * v ^ 2 := mul_le_mul_of_nonneg_right (le_abs_self y) (sq_nonneg v)
              _ ≤ |y| * 1 := mul_le_mul_of_nonneg_left hv2 (abs_nonneg y)
              _ ≤ |x| + 1 := by linarith
        _ = exp (|x| + 1) := one_mul _))
    (integrable_const _)
    (Eventually.of_forall fun v y _ => by
      have h : HasDerivAt (fun y : ℝ => exp (y * v ^ 2)) (exp (y * v ^ 2) * v ^ 2) y := by
        simpa using ((hasDerivAt_id y).mul_const (v ^ 2)).exp
      exact (h.const_mul (v ^ (2 * j))).congr_deriv (by ring))
  have hK : ∀ i : ℕ,
      kummerMoment i = fun y => ∫ v in Ioc (0 : ℝ) 1, v ^ (2 * i) * exp (y * v ^ 2) := by
    intro i; funext y; exact intervalIntegral.integral_of_le zero_le_one
  rw [hK j, hK (j + 1)]
  exact key.2

theorem kummerMoment_ode (j : ℕ) (x : ℝ) :
    2 * x * kummerMoment (j + 1) x + (2 * j + 1) * kummerMoment j x = exp x := by
  have hd : ∀ v ∈ uIcc (0 : ℝ) 1, HasDerivAt (fun v => v ^ (2 * j + 1) * exp (x * v ^ 2))
      ((2 * j + 1) * (v ^ (2 * j) * exp (x * v ^ 2)) +
        2 * x * (v ^ (2 * (j + 1)) * exp (x * v ^ 2))) v := by
    intro v _
    have h1 : HasDerivAt (fun v : ℝ => v ^ (2 * j + 1)) ((2 * j + 1 : ℝ) * v ^ (2 * j)) v := by
      simpa using hasDerivAt_pow (2 * j + 1) v
    have h2 : HasDerivAt (fun v : ℝ => exp (x * v ^ 2)) (exp (x * v ^ 2) * (x * (2 * v))) v := by
      simpa using ((hasDerivAt_pow 2 v).const_mul x).exp
    exact (h1.mul h2).congr_deriv (by ring)
  have hi : ∀ i : ℕ, IntervalIntegrable (fun v : ℝ => v ^ (2 * i) * exp (x * v ^ 2)) volume 0 1 :=
    fun i => (kummerMoment_continuous_integrand i x).intervalIntegrable _ _
  have := intervalIntegral.integral_eq_sub_of_hasDerivAt hd
    (((hi j).const_mul _).add ((hi (j + 1)).const_mul _))
  rw [intervalIntegral.integral_add ((hi j).const_mul _) ((hi (j + 1)).const_mul _),
    intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul] at this
  simp only [one_pow, mul_one, one_mul, ne_eq, add_eq_zero, one_ne_zero, and_false,
    not_false_eq_true, zero_pow, zero_mul, sub_zero] at this
  unfold kummerMoment
  linarith

/-- `₁F₁(1; 3/2; −x) = D(√x)/√x = e^{−x} K₀(x)`: Dawson's function in Kummer's form. -/
noncomputable def kummerHalf (x : ℝ) : ℝ := exp (-x) * kummerMoment 0 x

theorem kummerHalf_pos (x : ℝ) : 0 < kummerHalf x := mul_pos (exp_pos _) (kummerMoment_pos 0 x)

/-- `K_j(x) = Σ xᵏ/(k! (2k + 2j + 1))`, the series `ResidualShrinkage.LogKummer` and the checks
sum. -/
theorem kummerMoment_hasSum (j : ℕ) (x : ℝ) :
    HasSum (fun k : ℕ => x ^ k / (k.factorial * (2 * k + 2 * j + 1))) (kummerMoment j x) := by
  have hF : ∀ k : ℕ, (∫ v in (0 : ℝ)..1, v ^ (2 * j) * ((x * v ^ 2) ^ k / k.factorial)) =
      x ^ k / (k.factorial * (2 * k + 2 * j + 1)) := by
    intro k
    have he : (fun v : ℝ => v ^ (2 * j) * ((x * v ^ 2) ^ k / k.factorial)) =
        fun v => (x ^ k / k.factorial) * v ^ (2 * k + 2 * j) := by
      funext v; rw [mul_pow, ← pow_mul]; ring
    rw [he, intervalIntegral.integral_const_mul, integral_pow]
    have : (0 : ℝ) < k.factorial := by exact_mod_cast k.factorial_pos
    push_cast
    field_simp
    ring
  have hs := intervalIntegral.hasSum_integral_of_dominated_convergence (μ := volume)
    (a := 0) (b := 1)
    (F := fun k v => v ^ (2 * j) * ((x * v ^ 2) ^ k / k.factorial))
    (f := fun v => v ^ (2 * j) * exp (x * v ^ 2))
    (fun (k : ℕ) (_ : ℝ) => |x| ^ k / (k.factorial : ℝ))
    (fun k => (by fun_prop : Continuous fun v : ℝ =>
      v ^ (2 * j) * ((x * v ^ 2) ^ k / k.factorial)).aestronglyMeasurable)
    (fun k => Eventually.of_forall fun v hv => by
      rw [Set.uIoc_of_le zero_le_one] at hv
      have h1 : |v| ≤ 1 := by rw [abs_of_pos hv.1]; exact hv.2
      rw [Real.norm_eq_abs, abs_mul, abs_div, abs_pow, abs_pow, abs_mul, abs_pow,
        abs_of_pos (by exact_mod_cast k.factorial_pos : (0 : ℝ) < k.factorial)]
      have h2 : |v| ^ 2 ≤ 1 := pow_le_one₀ (abs_nonneg v) h1
      have h3 : (|x| * |v| ^ 2) ^ k ≤ |x| ^ k := by
        rw [mul_pow]
        exact mul_le_of_le_one_right (pow_nonneg (abs_nonneg x) _)
          (pow_le_one₀ (pow_nonneg (abs_nonneg v) _) h2)
      have hf : (0 : ℝ) < k.factorial := by exact_mod_cast k.factorial_pos
      calc |v| ^ (2 * j) * ((|x| * |v| ^ 2) ^ k / k.factorial) ≤ 1 * (|x| ^ k / k.factorial) :=
            mul_le_mul (pow_le_one₀ (abs_nonneg v) h1)
              (div_le_div_of_nonneg_right h3 hf.le) (by positivity) zero_le_one
        _ = |x| ^ k / k.factorial := one_mul _)
    (Eventually.of_forall fun _ _ => by simpa using Real.summable_pow_div_factorial |x|)
    intervalIntegrable_const
    (Eventually.of_forall fun v _ => by
      have := NormedSpace.expSeries_div_hasSum_exp (𝔸 := ℝ) (x * v ^ 2)
      rw [← Real.exp_eq_exp_ℝ] at this
      exact this.mul_left (v ^ (2 * j)))
  simp only [hF] at hs
  exact hs

/-- `g(x) = K₁(x)/K₀(x)`, the tilted mean of `V²`. -/
noncomputable def robasG (x : ℝ) : ℝ := kummerMoment 1 x / kummerMoment 0 x

theorem robasG_pos (x : ℝ) : 0 < robasG x := div_pos (kummerMoment_pos 1 x) (kummerMoment_pos 0 x)

theorem robasG_le_one (x : ℝ) : robasG x ≤ 1 :=
  (div_le_one (kummerMoment_pos 0 x)).2 (kummerMoment_succ_le 0 x)

theorem hasDerivAt_robasG (x : ℝ) :
    HasDerivAt robasG (kummerMoment 2 x / kummerMoment 0 x - robasG x ^ 2) x := by
  have h0 := (kummerMoment_pos 0 x).ne'
  exact ((hasDerivAt_kummerMoment 1 x).div (hasDerivAt_kummerMoment 0 x) h0).congr_deriv (by
    unfold robasG; field_simp)

/-- Cauchy–Schwarz for the moments: `K₁² ≤ K₀ K₂`. -/
theorem kummerMoment_sq_le (x : ℝ) :
    kummerMoment 1 x ^ 2 ≤ kummerMoment 0 x * kummerMoment 2 x := by
  have hi : ∀ i : ℕ, IntervalIntegrable (fun v : ℝ => v ^ (2 * i) * exp (x * v ^ 2)) volume 0 1 :=
    fun i => (kummerMoment_continuous_integrand i x).intervalIntegrable _ _
  have hq : ∀ l : ℝ,
      0 ≤ l ^ 2 * kummerMoment 0 x + 2 * l * kummerMoment 1 x + kummerMoment 2 x := by
    intro l
    have h0 : 0 ≤ ∫ v in (0 : ℝ)..1, (l + v ^ 2) ^ 2 * exp (x * v ^ 2) :=
      intervalIntegral.integral_nonneg zero_le_one fun v _ =>
        mul_nonneg (sq_nonneg _) (exp_pos _).le
    have he : (fun v : ℝ => (l + v ^ 2) ^ 2 * exp (x * v ^ 2)) = fun v =>
        l ^ 2 * (v ^ (2 * 0) * exp (x * v ^ 2)) + 2 * l * (v ^ (2 * 1) * exp (x * v ^ 2)) +
          v ^ (2 * 2) * exp (x * v ^ 2) := by funext v; ring
    rw [he, intervalIntegral.integral_add (((hi 0).const_mul _).add ((hi 1).const_mul _)) (hi 2),
      intervalIntegral.integral_add ((hi 0).const_mul _) ((hi 1).const_mul _),
      intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul] at h0
    unfold kummerMoment
    linarith
  have hK := kummerMoment_pos 0 x
  have := hq (-(kummerMoment 1 x / kummerMoment 0 x))
  field_simp at this
  nlinarith [this]

theorem robasG_monotone : Monotone robasG := by
  refine monotone_of_deriv_nonneg (fun x => (hasDerivAt_robasG x).differentiableAt) fun x => ?_
  rw [(hasDerivAt_robasG x).deriv]
  have hK := kummerMoment_pos 0 x
  have hs := kummerMoment_sq_le x
  unfold robasG
  rw [div_pow, sub_nonneg, div_le_div_iff₀ (pow_pos hK 2) hK]
  nlinarith

/-- `Ω(y) = −2 ln ₁F₁(1; 3/2; −y²/2)`. -/
noncomputable def robasBigOmega (y : ℝ) : ℝ := -2 * log (kummerHalf (y ^ 2 / 2))

/-- `ω = Ω′ = 2y(1 − g(y²/2))`. -/
noncomputable def robasOmega (y : ℝ) : ℝ := 2 * y * (1 - robasG (y ^ 2 / 2))

/-- `ω′ = 4g(1 − x(1 − g))` at `x = y²/2`. -/
noncomputable def robasOmegaDeriv (y : ℝ) : ℝ :=
  4 * robasG (y ^ 2 / 2) * (1 - y ^ 2 / 2 * (1 - robasG (y ^ 2 / 2)))

/-- `c = ω − y ω′`. -/
noncomputable def robasC (y : ℝ) : ℝ := robasOmega y - y * robasOmegaDeriv y

/-- `κ = 2Ω − y ω`, with `κ′ = c`. -/
noncomputable def robasKappa (y : ℝ) : ℝ := 2 * robasBigOmega y - y * robasOmega y

theorem robasBigOmega_eq (y : ℝ) :
    robasBigOmega y = y ^ 2 - 2 * log (kummerMoment 0 (y ^ 2 / 2)) := by
  unfold robasBigOmega kummerHalf
  rw [log_mul (exp_pos _).ne' (kummerMoment_pos 0 _).ne', log_exp]
  ring

theorem hasDerivAt_half_sq (y : ℝ) : HasDerivAt (fun y : ℝ => y ^ 2 / 2) y y := by
  simpa using (hasDerivAt_pow 2 y).div_const 2

theorem hasDerivAt_robasBigOmega (y : ℝ) : HasDerivAt robasBigOmega (robasOmega y) y := by
  have hK := (hasDerivAt_kummerMoment 0 (y ^ 2 / 2)).comp y (hasDerivAt_half_sq y)
  have hL := hK.log (kummerMoment_pos 0 _).ne'
  rw [show robasBigOmega = fun y => y ^ 2 - 2 * log (kummerMoment 0 (y ^ 2 / 2)) from
    funext robasBigOmega_eq]
  refine ((hasDerivAt_pow 2 y).sub (hL.const_mul 2)).congr_deriv ?_
  have h0 := (kummerMoment_pos 0 (y ^ 2 / 2)).ne'
  unfold robasOmega robasG
  simp only [Function.comp, zero_add, Nat.cast_ofNat]
  field_simp
  ring

theorem hasDerivAt_robasOmega (y : ℝ) : HasDerivAt robasOmega (robasOmegaDeriv y) y := by
  have hg := (hasDerivAt_robasG (y ^ 2 / 2)).comp y (hasDerivAt_half_sq y)
  refine (((hasDerivAt_id y).const_mul 2).mul ((hasDerivAt_const y (1 : ℝ)).sub hg)).congr_deriv ?_
  have o0 := kummerMoment_ode 0 (y ^ 2 / 2)
  have o1 := kummerMoment_ode 1 (y ^ 2 / 2)
  have h0 := (kummerMoment_pos 0 (y ^ 2 / 2)).ne'
  simp only [zero_add, Nat.reduceAdd] at o0 o1
  simp only [Pi.sub_apply, Function.comp_apply, id]
  unfold robasOmegaDeriv robasG
  field_simp
  linear_combination (4 * kummerMoment 0 (y ^ 2 / 2)) * o0 - (4 * kummerMoment 0 (y ^ 2 / 2)) * o1

theorem hasDerivAt_robasKappa (y : ℝ) : HasDerivAt robasKappa (robasC y) y :=
  (((hasDerivAt_robasBigOmega y).const_mul 2).sub
    ((hasDerivAt_id y).mul (hasDerivAt_robasOmega y))).congr_deriv (by unfold robasC; simp; ring)

theorem robasOmegaDeriv_le (y : ℝ) : robasOmegaDeriv y ≤ 4 := by
  unfold robasOmegaDeriv
  have hg := robasG_pos (y ^ 2 / 2)
  have hg1 := robasG_le_one (y ^ 2 / 2)
  have hx : 0 ≤ y ^ 2 / 2 * (1 - robasG (y ^ 2 / 2)) := mul_nonneg (by positivity) (by linarith)
  nlinarith

/-! #### The certified constants and the four bounds

`conformal-checks.py` reads these constants from this file and certifies `RobasDawsonBounds` with
them. -/

/-- An upper bound on `|ω|`. -/
noncomputable def robasOmegaMax : ℝ := 18841 / 10000

/-- A lower bound on `ω′`. -/
noncomputable def robasOmegaDerivMin : ℝ := -5511 / 10000

/-- An upper bound on `|c|`. -/
noncomputable def robasCMax : ℝ := 32365 / 10000

/-- The long window: at least `√2`. -/
noncomputable def robasWindowLong : ℝ := 1415 / 1000

/-- An upper bound on `κ`'s increment over the long window. -/
noncomputable def robasWindowLongBound : ℝ := 42766 / 10000

/-- The short window. -/
noncomputable def robasWindowShort : ℝ := 1351 / 1000

/-- An upper bound on `κ`'s increment over the short window. -/
noncomputable def robasWindowShortBound : ℝ := 41060 / 10000

/-- The four bounds on the module's Dawson function that RoBAS-Full's interval theorem takes,
certified by interval arithmetic in `conformal-checks.py`: `|ω| ≤ 1.8841`, `ω′ ≥ −0.5511`,
`|c| ≤ 3.2365`, and `κ`'s increments over windows of 1.415 and 1.351 at most 4.2766 and 4.1060. -/
structure RobasDawsonBounds : Prop where
  omega : ∀ y, |robasOmega y| ≤ robasOmegaMax
  omegaDeriv : ∀ y, robasOmegaDerivMin ≤ robasOmegaDeriv y
  c : ∀ y, |robasC y| ≤ robasCMax
  windowLong : ∀ a b, |b - a| ≤ robasWindowLong →
    |robasKappa b - robasKappa a| ≤ robasWindowLongBound
  windowShort : ∀ a b, |b - a| ≤ robasWindowShort →
    |robasKappa b - robasKappa a| ≤ robasWindowShortBound

/-- The mean value bound: a derivative at most `M` in size moves `f` by at most `M |b − a|`. -/
theorem abs_sub_le_of_hasDerivAt {f f' : ℝ → ℝ} {M : ℝ} (hf : ∀ y, HasDerivAt f (f' y) y)
    (hM : ∀ y, |f' y| ≤ M) (a b : ℝ) : |f b - f a| ≤ M * |b - a| := by
  have key : ∀ a b : ℝ, a < b → |f b - f a| ≤ M * |b - a| := by
    intro a b h
    obtain ⟨ξ, -, hξ⟩ := exists_hasDerivAt_eq_slope f f' h
      (fun y _ => (hf y).continuousAt.continuousWithinAt) fun y _ => hf y
    rw [eq_div_iff (sub_ne_zero.2 h.ne')] at hξ
    rw [← hξ, abs_mul]
    exact mul_le_mul_of_nonneg_right (hM ξ) (abs_nonneg _)
  rcases lt_trichotomy a b with h | h | h
  · exact key a b h
  · subst h; simp
  · rw [abs_sub_comm, abs_sub_comm b a]; exact key b a h

/-- Taylor's lower bound from `ω′ ≥ k`: `Ω q ≥ Ω p + ω p (q − p) + k (q − p)²/2`. -/
theorem robasBigOmega_taylor {k : ℝ} (hk : ∀ y, k ≤ robasOmegaDeriv y) (p q : ℝ) :
    robasBigOmega p + robasOmega p * (q - p) + k / 2 * (q - p) ^ 2 ≤ robasBigOmega q := by
  set φ : ℝ → ℝ := fun y => robasBigOmega y - k / 2 * y ^ 2
  have hφ : ∀ y, HasDerivAt φ (robasOmega y - k * y) y := fun y =>
    ((hasDerivAt_robasBigOmega y).sub ((hasDerivAt_pow 2 y).const_mul (k / 2))).congr_deriv
      (by push_cast; ring)
  have hψ : ∀ y, HasDerivAt (fun y => robasOmega y - k * y) (robasOmegaDeriv y - k) y := fun y =>
    ((hasDerivAt_robasOmega y).sub ((hasDerivAt_id y).const_mul k)).congr_deriv (by simp)
  have hmono : Monotone fun y => robasOmega y - k * y :=
    monotone_of_deriv_nonneg (fun y => (hψ y).differentiableAt) fun y => by
      rw [(hψ y).deriv]; linarith [hk y]
  have key : ∀ a b : ℝ, a < b → ∃ ξ ∈ Set.Ioo a b,
      φ b - φ a = (robasOmega ξ - k * ξ) * (b - a) := by
    intro a b h
    obtain ⟨ξ, hξ, hs⟩ := exists_hasDerivAt_eq_slope φ (fun y => robasOmega y - k * y) h
      (fun y _ => (hφ y).continuousAt.continuousWithinAt) fun y _ => hφ y
    refine ⟨ξ, hξ, ?_⟩
    rw [hs]; field_simp [sub_ne_zero.2 h.ne']
  rcases lt_trichotomy p q with h | h | h
  · obtain ⟨ξ, hξ, hs⟩ := key p q h
    have hm := hmono hξ.1.le
    simp only [φ] at hs hm
    nlinarith
  · subst h; simp
  · obtain ⟨ξ, hξ, hs⟩ := key q p h
    have hm := hmono hξ.2.le
    simp only [φ] at hs hm
    nlinarith

theorem abs_le_of_sq_le_sq'' {a b : ℝ} (h : a ^ 2 ≤ b ^ 2) (hb : 0 ≤ b) : |a| ≤ b := by
  by_contra hc
  push Not at hc
  have : b ^ 2 < |a| ^ 2 := pow_lt_pow_left₀ hc hb (by norm_num)
  rw [sq_abs] at this
  linarith

/-- The slope at a zero, from the four bounds: with `d = p − q`, `t − u = d`,
`(n + 1)(t + u) d = −(Ω p − Ω q)` and the other units' deviations admissible,
`(n − 1) d² + (d ω p − (Ω p − Ω q)) + t d (κ p − κ q) > 0` for `n ≥ 4`. -/
theorem robas_core (hb : RobasDawsonBounds) {n : ℕ} (hn : 4 ≤ n) {p q t u : ℝ} (hpq : p ≠ q)
    (htu : t - u = p - q)
    (hsum : ((n : ℝ) + 1) * (t + u) * (p - q) = -(robasBigOmega p - robasBigOmega q))
    (hadm : (t + u) ^ 2 ≤ ((n : ℝ) - 1) * (1 - t ^ 2 - u ^ 2)) :
    0 < ((n : ℝ) - 1) * (p - q) ^ 2 +
      ((p - q) * robasOmega p - (robasBigOmega p - robasBigOmega q)) +
        t * (p - q) * (robasKappa p - robasKappa q) := by
  have hn' : (4 : ℝ) ≤ n := by exact_mod_cast hn
  set d := p - q with hd_def
  set E := robasBigOmega p - robasBigOmega q with hE_def
  set K := robasKappa p - robasKappa q with hK_def
  set s := t + u with hs_def
  have hd0 : d ≠ 0 := sub_ne_zero.2 hpq
  have hdpos : 0 < |d| := abs_pos.2 hd0
  have hd2 : 0 < d ^ 2 := by positivity
  have hm : (0 : ℝ) < n + 1 := by linarith
  -- Ω̄ = E/d = −(n + 1) s
  have hE : E = -((n + 1) * s) * d := by linarith [hsum]
  have htd : t * d = (d ^ 2 + s * d) / 2 := by
    have : t = (d + s) / 2 := by linarith [htu]
    rw [this]; ring
  -- the three pointwise bounds
  have h1 : robasOmegaDerivMin / 2 * d ^ 2 ≤ d * robasOmega p - E := by
    have hT := robasBigOmega_taylor hb.omegaDeriv p q
    have : q - p = -d := by rw [hd_def]; ring
    rw [this] at hT
    simp only [hE_def]
    nlinarith [hT]
  have hEb : |E| ≤ robasOmegaMax * |d| := by
    have := abs_sub_le_of_hasDerivAt hasDerivAt_robasBigOmega hb.omega q p
    rwa [← hE_def, ← hd_def] at this
  have hKc : |K| ≤ robasCMax * |d| := by
    have := abs_sub_le_of_hasDerivAt hasDerivAt_robasKappa hb.c q p
    rwa [← hK_def, ← hd_def] at this
  have hOb : (n + 1) * |s| ≤ robasOmegaMax := by
    have h' : |E| = (n + 1) * |s| * |d| := by
      rw [hE, abs_mul, abs_neg, abs_mul, abs_of_pos hm]
    rw [h'] at hEb
    exact le_of_mul_le_mul_right hEb hdpos
  -- admissibility: (n − 1)(s² + d²)/2 ≤ (n − 1) − s²
  have hadm' : s ^ 2 ≤ ((n : ℝ) - 1) * (1 - (s ^ 2 + d ^ 2) / 2) := by
    have : t ^ 2 + u ^ 2 = (s ^ 2 + d ^ 2) / 2 := by
      have ht : t = (d + s) / 2 := by linarith [htu]
      have hu : u = (s - d) / 2 := by linarith [htu]
      rw [ht, hu]; ring
    nlinarith [hadm, this]
  have hd2le : d ^ 2 ≤ 2 := by nlinarith [sq_nonneg s]
  have hdL : |d| ≤ robasWindowLong := by
    refine abs_le_of_sq_le_sq'' ?_ (by norm_num [robasWindowLong])
    norm_num [robasWindowLong]; linarith
  -- Ξ ≥ d² ((n − 1) + ω̲/2 − |K|/2 − (n + 1)|s| c̄/2)
  have hlow : d ^ 2 * (((n : ℝ) - 1) + robasOmegaDerivMin / 2 - |K| / 2 -
      |s| * robasCMax / 2) ≤
      ((n : ℝ) - 1) * d ^ 2 + (d * robasOmega p - E) + t * d * K := by
    rw [htd]
    have hA : -(|K|) ≤ K := neg_abs_le K
    have hB : |s * d * K| ≤ |s| * (robasCMax * |d| * |d|) := by
      rw [abs_mul, abs_mul]
      calc |s| * |d| * |K| = |s| * (|K| * |d|) := by ring
        _ ≤ |s| * (robasCMax * |d| * |d|) :=
          mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_right hKc (abs_nonneg _)) (abs_nonneg _)
    have hB' : -(|s| * (robasCMax * |d| * |d|)) ≤ s * d * K :=
      le_trans (neg_le_neg hB) (neg_abs_le _)
    have hdd : |d| * |d| = d ^ 2 := by rw [← abs_mul, ← sq, abs_of_nonneg (sq_nonneg d)]
    have hk2 : d ^ 2 / 2 * (-(|K|)) ≤ d ^ 2 / 2 * K := mul_le_mul_of_nonneg_left hA (by positivity)
    have hB'' : -(|s| * robasCMax * d ^ 2) ≤ s * d * K := by
      rw [← hdd]; linarith [hB']
    linarith [h1, hk2, hB'']
  refine lt_of_lt_of_le (mul_pos hd2 ?_) hlow
  simp only [robasOmegaDerivMin, robasCMax]
  by_cases h5 : 5 ≤ n
  · have hn5 : (5 : ℝ) ≤ n := by exact_mod_cast h5
    have hKW := hb.windowLong q p hdL
    have hs6 : 6 * |s| ≤ robasOmegaMax := by nlinarith only [hOb, hn5, abs_nonneg s]
    simp only [robasWindowLongBound, robasOmegaMax] at hKW hs6
    linarith only [hKW, hs6, hn5]
  · have hn4 : (n : ℝ) = 4 := by
      have : n = 4 := by omega
      rw [this]; norm_num
    rw [hn4] at hOb hadm'
    rw [hn4]
    by_cases hc : (5 * s) ^ 2 ≤ 15 * (2 - robasWindowShort ^ 2)
    · have hKW := hb.windowLong q p hdL
      have hs5 : 5 * |s| ≤ 16193 / 10000 := by
        have : |5 * s| ≤ 16193 / 10000 := abs_le_of_sq_le_sq'' (by
          simp only [robasWindowShort] at hc; norm_num at hc ⊢; linarith only [hc]) (by norm_num)
        rwa [abs_mul, abs_of_pos (by norm_num : (0 : ℝ) < 5)] at this
      simp only [robasWindowLongBound] at hKW
      linarith only [hKW, hs5]
    · push Not at hc
      have hdS : |d| ≤ robasWindowShort := by
        refine abs_le_of_sq_le_sq'' ?_ (by norm_num [robasWindowShort])
        simp only [robasWindowShort] at hc ⊢
        have h25 : (5 * s) ^ 2 = 25 * s ^ 2 := by ring
        linarith only [hadm', hc, h25]
      have hKW := hb.windowShort q p hdS
      have hs5 : 5 * |s| ≤ robasOmegaMax := by linarith only [hOb]
      simp only [robasWindowShortBound, robasOmegaMax] at hKW hs5
      linarith only [hKW, hs5]

/-! #### RoBAS-Full's score, its p-value, and the score difference -/

/-- RoBAS-Full's score of a unit whose `n` others have sum `so` and sum of squares `qo`, at the
plug-in `σ̂² = s2`: the log marginal likelihood `−n s²/(2σ̂²) + ln ₁F₁(1; 3/2; −n r̄²/(2σ̂²))` of
the others under the horseshoe residual model, up to the constant every unit shares (RoBAS's
RoBAS-Full derivation; `ResidualShrinkage.Score`'s Horseshoe branch). Here `r̄ = so/n` and
`s² = qo/n − r̄²`. -/
noncomputable def robasFullScore (n : ℕ) (s2 so qo : ℝ) : ℝ :=
  -(n : ℝ) * (qo / n - (so / n) ^ 2) / (2 * s2) +
    log (kummerHalf ((n : ℝ) * (so / n) ^ 2 / (2 * s2)))

theorem two_mul_robasFullScore {n : ℕ} (hn : (n : ℝ) ≠ 0) {T : ℝ} (hT : 0 < T) (so qo : ℝ) :
    2 * T * robasFullScore n (T / n) so qo =
      -(n : ℝ) * qo + so ^ 2 - T * robasBigOmega (so / sqrt T) := by
  unfold robasFullScore robasBigOmega
  have hs : sqrt T ^ 2 = T := sq_sqrt hT.le
  have hx : (n : ℝ) * (so / n) ^ 2 / (2 * (T / n)) = (so / sqrt T) ^ 2 / 2 := by
    rw [div_pow so (sqrt T), hs]; field_simp
  rw [hx]
  field_simp
  ring

variable {n : ℕ}

/-- All `n + 1` residuals' sum of squared deviations, `n σ̂²`. -/
noncomputable def robasSpread (r : Fin n → ℝ) (z : ℝ) : ℝ :=
  (∑ i, r i ^ 2 + z ^ 2) - (∑ i, r i + z) ^ 2 / (n + 1)

/-- Calibration unit `j`'s RoBAS-Full score with candidate residual `z`: its others are the other
calibration residuals and `z`. -/
noncomputable def robasUnitScore (r : Fin n → ℝ) (z : ℝ) (j : Fin n) : ℝ :=
  robasFullScore n (robasSpread r z / n) (∑ i, r i + z - r j) (∑ i, r i ^ 2 + z ^ 2 - r j ^ 2)

/-- The candidate's RoBAS-Full score: its others are the calibration residuals. -/
noncomputable def robasCandidateScore (r : Fin n → ℝ) (z : ℝ) : ℝ :=
  robasFullScore n (robasSpread r z / n) (∑ i, r i) (∑ i, r i ^ 2)

/-- The conformal p-value of candidate residual `z`: the share of the `n + 1` units whose score is
at least the candidate's, its own tie counted (`ResidualShrinkage`'s `atOrAbove + 1`). -/
noncomputable def robasPValue (r : Fin n → ℝ) (z : ℝ) : ℝ :=
  (((Finset.univ.filter fun j => robasCandidateScore r z ≤ robasUnitScore r z j).card : ℝ) + 1) /
    (n + 1)

/-- `Δ_j(z) = (r_j − z) L_j(z) − T(z) (Ω(p) − Ω(q))`, with `L_j(z) = (n − 1)z + (n + 1)r_j − 2S` and
`p`, `q` the two leave-one-out sums over `√T`. -/
noncomputable def robasDelta (r : Fin n → ℝ) (j : Fin n) (z : ℝ) : ℝ :=
  (r j - z) * ((n - 1) * z + (n + 1) * r j - 2 * ∑ i, r i) -
    robasSpread r z * (robasBigOmega ((∑ i, r i + z - r j) / sqrt (robasSpread r z)) -
      robasBigOmega ((∑ i, r i) / sqrt (robasSpread r z)))

theorem robas_score_diff (hn : (n : ℝ) ≠ 0) (r : Fin n → ℝ) (z : ℝ) (hT : 0 < robasSpread r z)
    (j : Fin n) :
    2 * robasSpread r z * (robasUnitScore r z j - robasCandidateScore r z) = robasDelta r j z := by
  unfold robasUnitScore robasCandidateScore robasDelta
  rw [mul_sub, two_mul_robasFullScore hn hT, two_mul_robasFullScore hn hT]
  ring

theorem robasSpread_eq (r : Fin n → ℝ) (z : ℝ) : robasSpread r z =
    ∑ i, (r i - (∑ k, r k + z) / (n + 1)) ^ 2 + (z - (∑ k, r k + z) / (n + 1)) ^ 2 := by
  unfold robasSpread
  set μ := (∑ k, r k + z) / (n + 1) with hμ
  have hm : ((n : ℝ) + 1) ≠ 0 := by positivity
  have hsum : ∑ i, (r i - μ) ^ 2 = ∑ i, r i ^ 2 - 2 * μ * ∑ i, r i + n * μ ^ 2 := by
    simp only [sub_sq, Finset.sum_add_distrib, Finset.sum_sub_distrib, Finset.sum_const,
      Finset.card_univ, Fintype.card_fin, nsmul_eq_mul, ← Finset.mul_sum, ← Finset.sum_mul]
    ring
  rw [hsum, hμ]
  field_simp
  ring

theorem robasSpread_pos (r : Fin n → ℝ) (hr : ∃ i k, r i ≠ r k) (z : ℝ) : 0 < robasSpread r z := by
  obtain ⟨i, k, hik⟩ := hr
  have hne : i ≠ k := fun h => hik (h ▸ rfl)
  rw [robasSpread_eq]
  set μ := (∑ k, r k + z) / (n + 1)
  have h2 : (r i - μ) ^ 2 + (r k - μ) ^ 2 ≤ ∑ l, (r l - μ) ^ 2 := by
    have := Finset.sum_le_sum_of_subset_of_nonneg (Finset.subset_univ {i, k})
      (fun l _ _ => sq_nonneg (r l - μ))
    rwa [Finset.sum_pair hne] at this
  have h3 : 0 < (r i - μ) ^ 2 + (r k - μ) ^ 2 := by
    by_contra h
    push Not at h
    have a1 : (r i - μ) ^ 2 = 0 := le_antisymm (by nlinarith [sq_nonneg (r k - μ)]) (sq_nonneg _)
    have a2 : (r k - μ) ^ 2 = 0 := le_antisymm (by nlinarith [sq_nonneg (r i - μ)]) (sq_nonneg _)
    have := pow_eq_zero_iff (two_ne_zero) |>.1 a1
    have := pow_eq_zero_iff (two_ne_zero) |>.1 a2
    exact hik (by linarith)
  linarith [sq_nonneg (z - μ)]

/-- The slope at a zero, as algebra: with `σ² = T`, `p`, `q` the leave-one-out sums over `σ`,
`t` the candidate's deviation over `σ`, and `Δ = 0`, the derivative times `a − z` is `T` times
the core quantity. -/
theorem robas_slope_algebra {m S a z σ T E ωp ωq : ℝ} (hσ : 0 < σ) (hT : σ ^ 2 = T)
    (hm : 0 < m)
    (hz : (a - z) * ((m - 2) * z + m * a - 2 * S) - T * E = 0) :
    (-((m - 2) * z + m * a - 2 * S) + (a - z) * (m - 2) -
        ((2 * z - 2 * (S + z) / m) * E +
          T * (ωp * ((1 * σ - (S + z - a) * ((2 * z - 2 * (S + z) / m) / (2 * σ))) / σ ^ 2) -
            ωq * ((0 * σ - S * ((2 * z - 2 * (S + z) / m) / (2 * σ))) / σ ^ 2)))) * (a - z) =
      T * ((m - 2) * ((S + z - a) / σ - S / σ) ^ 2 +
        (((S + z - a) / σ - S / σ) * ωp - E) +
          (z - (S + z) / m) / σ * ((S + z - a) / σ - S / σ) *
            (2 * E - ((S + z - a) / σ * ωp - S / σ * ωq))) := by
  subst hT
  field_simp
  linear_combination (-(m * σ)) * hz

theorem hasDerivAt_robasSpread (r : Fin n → ℝ) (z : ℝ) :
    HasDerivAt (robasSpread r) (2 * z - 2 * (∑ i, r i + z) / (n + 1)) z := by
  have h := ((hasDerivAt_pow 2 z).const_add (∑ i, r i ^ 2)).sub
    ((((hasDerivAt_id z).const_add (∑ i, r i)).pow 2).div_const ((n : ℝ) + 1))
  refine h.congr_deriv ?_
  simp only [id, Nat.cast_ofNat]
  ring

/-- `Δ_j`'s derivative. -/
theorem hasDerivAt_robasDelta (r : Fin n → ℝ) (hr : ∃ i k, r i ≠ r k) (j : Fin n) (z : ℝ) :
    HasDerivAt (robasDelta r j)
      ((0 - 1) * (((n : ℝ) - 1) * z + (n + 1) * r j - 2 * ∑ i, r i) + (r j - z) * ((n - 1) * 1) -
        ((2 * z - 2 * (∑ i, r i + z) / (n + 1)) *
            (robasBigOmega ((∑ i, r i + z - r j) / sqrt (robasSpread r z)) -
              robasBigOmega ((∑ i, r i) / sqrt (robasSpread r z))) +
          robasSpread r z *
            (robasOmega ((∑ i, r i + z - r j) / sqrt (robasSpread r z)) *
                ((1 * sqrt (robasSpread r z) - (∑ i, r i + z - r j) *
                    ((2 * z - 2 * (∑ i, r i + z) / (n + 1)) / (2 * sqrt (robasSpread r z)))) /
                  sqrt (robasSpread r z) ^ 2) -
              robasOmega ((∑ i, r i) / sqrt (robasSpread r z)) *
                ((0 * sqrt (robasSpread r z) - (∑ i, r i) *
                    ((2 * z - 2 * (∑ i, r i + z) / (n + 1)) / (2 * sqrt (robasSpread r z)))) /
                  sqrt (robasSpread r z) ^ 2)))) z := by
  have hT := robasSpread_pos r hr z
  have hσd := (hasDerivAt_robasSpread r z).sqrt hT.ne'
  have hσ0 : sqrt (robasSpread r z) ≠ 0 := (sqrt_pos.2 hT).ne'
  have hPd := ((((hasDerivAt_id z).const_add (∑ i, r i)).sub_const (r j)).div hσd hσ0)
  have hQd := ((hasDerivAt_const z (∑ i, r i)).div hσd hσ0)
  have hΩP := (hasDerivAt_robasBigOmega _).comp z hPd
  have hΩQ := (hasDerivAt_robasBigOmega _).comp z hQd
  have hL := (((hasDerivAt_id z).const_mul ((n : ℝ) - 1)).add_const ((n + 1) * r j)).sub_const
    (2 * ∑ i, r i)
  exact (((hasDerivAt_const z (r j)).sub (hasDerivAt_id z)).mul hL).sub
    ((hasDerivAt_robasSpread r z).mul (hΩP.sub hΩQ))

theorem robasDelta_self (r : Fin n → ℝ) (j : Fin n) : robasDelta r j (r j) = 0 := by
  unfold robasDelta
  simp

/-- Unit `j`'s `Δ_j` rises through each of its zeros left of `r_j` and falls through each right
of it. -/
theorem robasDelta_rises (hb : RobasDawsonBounds) (hn : 4 ≤ n) (r : Fin n → ℝ)
    (hr : ∃ i k, r i ≠ r k) (j : Fin n) (z : ℝ) (hz : robasDelta r j z = 0) (hzj : z ≠ r j) :
    ∃ D, HasDerivAt (robasDelta r j) D z ∧ 0 < D * (r j - z) := by
  refine ⟨_, hasDerivAt_robasDelta r hr j z, ?_⟩
  set S := ∑ i, r i with hS
  set a := r j with ha
  set T := robasSpread r z with hT_def
  have hT := robasSpread_pos r hr z
  set σ := sqrt T with hσ_def
  have hσ : 0 < σ := sqrt_pos.2 hT
  have hσ2 : σ ^ 2 = T := sq_sqrt hT.le
  have hm : (0 : ℝ) < n + 1 := by positivity
  have hn1 : (1 : ℝ) ≤ n := by exact_mod_cast (by omega : 1 ≤ n)
  set p := (S + z - a) / σ with hp
  set q := S / σ with hq
  set μ := (S + z) / (n + 1) with hμ
  set t := (z - μ) / σ with ht
  set u := (a - μ) / σ with hu
  have hz' : (a - z) * (((n : ℝ) + 1 - 2) * z + (n + 1) * a - 2 * S) -
      T * (robasBigOmega p - robasBigOmega q) = 0 := by
    rw [← hz]; unfold robasDelta; rw [← hS, ← ha, ← hT_def, ← hσ_def]; ring
  have halg := robas_slope_algebra (m := (n : ℝ) + 1) (S := S) (a := a) (z := z)
    (E := robasBigOmega p - robasBigOmega q) (ωp := robasOmega p) (ωq := robasOmega q) hσ hσ2 hm hz'
  -- the core's hypotheses
  have hpq : p ≠ q := by
    intro h
    apply hzj
    have : (S + z - a) / σ = S / σ := h
    field_simp at this
    linarith
  have htu : t - u = p - q := by simp only [ht, hu, hp, hq]; field_simp; ring
  have hsum : ((n : ℝ) + 1) * (t + u) * (p - q) = -(robasBigOmega p - robasBigOmega q) := by
    have hE : robasBigOmega p - robasBigOmega q =
        (a - z) * (((n : ℝ) + 1 - 2) * z + (n + 1) * a - 2 * S) / T := by
      have hT0 : T ≠ 0 := hT.ne'
      rw [eq_div_iff hT0]; linarith [hz']
    rw [hE, ← hσ2]
    simp only [ht, hu, hp, hq, hμ]
    field_simp
    ring
  have hadm : (t + u) ^ 2 ≤ ((n : ℝ) - 1) * (1 - t ^ 2 - u ^ 2) := by
    set e := Finset.univ.erase j
    have hcard : (e.card : ℝ) = n - 1 := by
      rw [Finset.card_erase_of_mem (Finset.mem_univ j), Finset.card_univ, Fintype.card_fin,
        Nat.cast_sub (by omega)]
      simp
    have hdev : ∑ i ∈ e, (r i - μ) = -((z - μ) + (a - μ)) := by
      have h1 := Finset.add_sum_erase Finset.univ (fun i => r i - μ) (Finset.mem_univ j)
      have h2 : ∑ i, (r i - μ) + (z - μ) = 0 := by
        rw [Finset.sum_sub_distrib, Finset.sum_const, Finset.card_univ, Fintype.card_fin,
          nsmul_eq_mul, hμ, ← hS]
        field_simp
        ring
      linarith
    have hsq : ∑ i ∈ e, (r i - μ) ^ 2 = T - (z - μ) ^ 2 - (a - μ) ^ 2 := by
      have h1 := Finset.add_sum_erase Finset.univ (fun i => (r i - μ) ^ 2) (Finset.mem_univ j)
      have h2 := robasSpread_eq r z
      rw [← hT_def, ← hS, ← hμ] at h2
      linarith
    have hcs := Finset.sum_mul_sq_le_sq_mul_sq e (fun i => r i - μ) (fun _ => (1 : ℝ))
    simp only [mul_one, one_pow, Finset.sum_const, nsmul_eq_mul] at hcs
    rw [hdev, hsq, hcard] at hcs
    simp only [ht, hu]
    rw [← add_div, div_pow, div_pow, div_pow, hσ2]
    have hT0 : T ≠ 0 := hT.ne'
    rw [div_le_iff₀ hT]
    have : ((n : ℝ) - 1) * (1 - (z - μ) ^ 2 / T - (a - μ) ^ 2 / T) * T =
        ((n : ℝ) - 1) * (T - (z - μ) ^ 2 - (a - μ) ^ 2) := by field_simp
    rw [this]
    calc (z - μ + (a - μ)) ^ 2 = (-(z - μ + (a - μ))) ^ 2 := by ring
      _ ≤ (T - (z - μ) ^ 2 - (a - μ) ^ 2) * ((n : ℝ) - 1) := hcs
      _ = ((n : ℝ) - 1) * (T - (z - μ) ^ 2 - (a - μ) ^ 2) := by ring
  have hcore := robas_core hb hn hpq htu hsum hadm
  have hκ : robasKappa p - robasKappa q =
      2 * (robasBigOmega p - robasBigOmega q) - (p * robasOmega p - q * robasOmega q) := by
    unfold robasKappa; ring
  have hfinal : 0 < T * (((n : ℝ) + 1 - 2) * (p - q) ^ 2 +
      ((p - q) * robasOmega p - (robasBigOmega p - robasBigOmega q)) +
        t * (p - q) * (2 * (robasBigOmega p - robasBigOmega q) -
          (p * robasOmega p - q * robasOmega q))) := by
    refine mul_pos hT ?_
    rw [← hκ]
    have : ((n : ℝ) + 1 - 2) = (n : ℝ) - 1 := by ring
    rw [this]
    exact hcore
  have hD : ((0 - 1) * (((n : ℝ) - 1) * z + (n + 1) * a - 2 * S) + (a - z) * ((n - 1) * 1) -
      ((2 * z - 2 * (S + z) / (n + 1)) * (robasBigOmega p - robasBigOmega q) +
        T * (robasOmega p * ((1 * σ - (S + z - a) * ((2 * z - 2 * (S + z) / (n + 1)) / (2 * σ))) /
            σ ^ 2) -
          robasOmega q * ((0 * σ - S * ((2 * z - 2 * (S + z) / (n + 1)) / (2 * σ))) / σ ^ 2)))) *
        (a - z) = T * (((n : ℝ) + 1 - 2) * (p - q) ^ 2 +
      ((p - q) * robasOmega p - (robasBigOmega p - robasBigOmega q)) +
        t * (p - q) * (2 * (robasBigOmega p - robasBigOmega q) -
          (p * robasOmega p - q * robasOmega q))) := by
    rw [← halg]
    ring
  rw [hD]
  exact hfinal

/-- A continuous function that rises through each of its zeros in `[a, b]` stays nonnegative from
`a`. -/
theorem nonneg_of_rising {f : ℝ → ℝ} {a b : ℝ} (hab : a ≤ b) (hf : ContinuousOn f (Set.Icc a b))
    (hz : ∀ x ∈ Set.Icc a b, f x = 0 → ∃ D, HasDerivAt f D x ∧ 0 < D) (ha : 0 ≤ f a) :
    0 ≤ f b := by
  by_contra hb
  push Not at hb
  set A := Set.Icc a b ∩ f ⁻¹' Set.Ici 0 with hA
  have hAc : IsClosed A := hf.preimage_isClosed_of_isClosed isClosed_Icc isClosed_Ici
  have hAne : A.Nonempty := ⟨a, Set.left_mem_Icc.2 hab, ha⟩
  have hAbdd : BddAbove A := ⟨b, fun x hx => hx.1.2⟩
  set s := sSup A
  have hsA : s ∈ A := hAc.csSup_mem hAne hAbdd
  have hsA2 : 0 ≤ f s := hsA.2
  have hsb : s < b := lt_of_le_of_ne hsA.1.2 fun h => by rw [h] at hsA2; linarith
  have hneg : ∀ x ∈ Set.Ioc s b, f x < 0 := fun x hx => by
    by_contra h
    push Not at h
    have : x ≤ s := le_csSup hAbdd ⟨⟨hsA.1.1.trans hx.1.le, hx.2⟩, h⟩
    linarith [hx.1]
  have hev : ∀ᶠ x in 𝓝[>] s, 0 < f x := by
    rcases eq_or_lt_of_le hsA2 with h0 | hpos
    · obtain ⟨D, hD, hDpos⟩ := hz s hsA.1 h0.symm
      have h2 : Tendsto (slope f s) (𝓝[>] s) (𝓝 D) :=
        (hasDerivAt_iff_tendsto_slope.1 hD).mono_left (nhdsWithin_mono _ fun x hx => ne_of_gt hx)
      filter_upwards [h2.eventually (lt_mem_nhds hDpos), self_mem_nhdsWithin] with x hx hxs
      rw [slope_def_field, ← h0, sub_zero] at hx
      exact (div_pos_iff_of_pos_right (sub_pos.2 hxs)).1 hx
    · have h3 : ∀ᶠ x in 𝓝[Set.Icc a b] s, 0 < f x := (hf s hsA.1).eventually (lt_mem_nhds hpos)
      have h4 : 𝓝[>] s ≤ 𝓝[Set.Icc a b] s := by
        rw [← nhdsWithin_Ioo_eq_nhdsGT hsb]
        exact nhdsWithin_mono _ fun x hx => ⟨hsA.1.1.trans hx.1.le, hx.2.le⟩
      exact h4 h3
  obtain ⟨x, hx1, hx2⟩ := (hev.and (Ioo_mem_nhdsGT hsb)).exists
  linarith [hneg x ⟨hx2.1, hx2.2.le⟩]

/-- A continuous function that falls through each of its zeros in `[a, b]` is nonnegative at `a`
when it is at `b`. -/
theorem nonneg_of_falling {f : ℝ → ℝ} {a b : ℝ} (hab : a ≤ b) (hf : ContinuousOn f (Set.Icc a b))
    (hz : ∀ x ∈ Set.Icc a b, f x = 0 → ∃ D, HasDerivAt f D x ∧ D < 0) (hb : 0 ≤ f b) :
    0 ≤ f a := by
  have key := nonneg_of_rising (f := fun x => f (-x)) (a := -b) (b := -a) (by linarith)
    (hf.comp continuousOn_neg fun x hx => ⟨by linarith [hx.2], by linarith [hx.1]⟩)
    (fun x hx h0 => by
      obtain ⟨D, hD, hDneg⟩ := hz (-x) ⟨by linarith [hx.2], by linarith [hx.1]⟩ h0
      exact ⟨D * -1, hD.comp x (hasDerivAt_neg x), by linarith⟩)
    (by simpa using hb)
  simpa using key

/-- Each unit's acceptance set `{z : Δ_j(z) ≥ 0}` is an interval. -/
theorem robas_acceptance_ordConnected (hb : RobasDawsonBounds) (hn : 4 ≤ n) (r : Fin n → ℝ)
    (hr : ∃ i k, r i ≠ r k) (j : Fin n) : Set.OrdConnected {z | 0 ≤ robasDelta r j z} := by
  have hcont : Continuous (robasDelta r j) :=
    continuous_iff_continuousAt.2 fun z => (hasDerivAt_robasDelta r hr j z).continuousAt
  refine ⟨fun x hx y hy w hw => ?_⟩
  simp only [Set.mem_setOf_eq] at hx hy ⊢
  rcases lt_trichotomy w (r j) with h | h | h
  · refine nonneg_of_rising hw.1 hcont.continuousOn (fun v hv h0 => ?_) hx
    have hv' : v < r j := lt_of_le_of_lt hv.2 h
    obtain ⟨D, hD, hpos⟩ := robasDelta_rises hb hn r hr j v h0 hv'.ne
    exact ⟨D, hD, pos_of_mul_pos_left hpos (sub_pos.2 hv').le⟩
  · rw [h, robasDelta_self]
  · refine nonneg_of_falling hw.2 hcont.continuousOn (fun v hv h0 => ?_) hy
    have hv' : r j < v := lt_of_lt_of_le h hv.1
    obtain ⟨D, hD, hpos⟩ := robasDelta_rises hb hn r hr j v h0 hv'.ne'
    exact ⟨D, hD, neg_of_mul_pos_left hpos (sub_neg.2 hv').le⟩

/-- A critical point of a function whose derivative is monotone is its minimum. -/
theorem le_of_monotone_deriv {F F' : ℝ → ℝ} (hF : ∀ x, HasDerivAt F (F' x) x) (hmono : Monotone F')
    {x₀ : ℝ} (h0 : F' x₀ = 0) (y : ℝ) : F x₀ ≤ F y := by
  rcases lt_trichotomy x₀ y with h | h | h
  · obtain ⟨ξ, hξ, hs⟩ := exists_hasDerivAt_eq_slope F F' h
      (fun x _ => (hF x).continuousAt.continuousWithinAt) fun x _ => hF x
    have := hmono hξ.1.le
    rw [h0, hs, le_div_iff₀ (sub_pos.2 h)] at this
    linarith
  · rw [h]
  · obtain ⟨ξ, hξ, hs⟩ := exists_hasDerivAt_eq_slope F F' h
      (fun x _ => (hF x).continuousAt.continuousWithinAt) fun x _ => hF x
    have := hmono hξ.2.le
    rw [h0, hs, div_le_iff₀ (sub_pos.2 h)] at this
    linarith

/-- A candidate residual at which every calibration unit's score is at least the candidate's: where
the candidate sits at the minimum of the scores' common convex profile. -/
theorem robas_common_point (r : Fin n → ℝ) (hr : ∃ i k, r i ≠ r k) :
    ∃ z₀, ∀ j, 0 ≤ robasDelta r j z₀ := by
  have hn : 1 ≤ n := by
    obtain ⟨i, -, -⟩ := hr
    exact Fin.pos_iff_nonempty.2 ⟨i⟩
  have hnR : (1 : ℝ) ≤ n := by exact_mod_cast hn
  set S := ∑ i, r i with hS
  have hT := robasSpread_pos r hr
  set h : ℝ → ℝ := fun z => n * z - S * robasG (S ^ 2 / (2 * robasSpread r z)) with hh
  have hTc : Continuous (robasSpread r) :=
    continuous_iff_continuousAt.2 fun z => (hasDerivAt_robasSpread r z).continuousAt
  have hgc : Continuous robasG :=
    continuous_iff_continuousAt.2 fun x => (hasDerivAt_robasG x).continuousAt
  have hc : Continuous h := by
    refine (continuous_const.mul continuous_id).sub (continuous_const.mul (hgc.comp ?_))
    exact continuous_const.div (continuous_const.mul hTc) fun z => by linarith [hT z]
  have hbound : ∀ z, |S * robasG (S ^ 2 / (2 * robasSpread r z))| ≤ |S| := fun z => by
    rw [abs_mul, abs_of_pos (robasG_pos _)]
    exact mul_le_of_le_one_right (abs_nonneg S) (robasG_le_one _)
  set M := |S| / n + 1
  have hM : 0 ≤ M := by positivity
  have hpos : 0 ≤ h M := by
    have := hbound M
    simp only [hh]
    have e : (n : ℝ) * M = |S| + n := by simp only [M]; field_simp
    rw [e]; linarith [abs_le.1 this]
  have hneg : h (-M) ≤ 0 := by
    have := hbound (-M)
    simp only [hh]
    have e : (n : ℝ) * -M = -(|S| + n) := by simp only [M]; field_simp
    rw [e]; linarith [abs_le.1 this]
  obtain ⟨z₀, -, hz₀⟩ := intermediate_value_Icc (by linarith) hc.continuousOn ⟨hneg, hpos⟩
  refine ⟨z₀, fun j => ?_⟩
  set T := robasSpread r z₀ with hT_def
  set σ := sqrt T with hσ_def
  have hσ : 0 < σ := sqrt_pos.2 (hT z₀)
  have hσ2 : σ ^ 2 = T := sq_sqrt (hT z₀).le
  set F : ℝ → ℝ := fun ρ => n * ρ ^ 2 + (S + z₀ - ρ) ^ 2 - T * robasBigOmega ((S + z₀ - ρ) / σ)
  set F' : ℝ → ℝ := fun ρ => 2 * n * ρ - 2 * (S + z₀ - ρ) + σ * robasOmega ((S + z₀ - ρ) / σ)
  have hin : ∀ ρ, HasDerivAt (fun ρ => (S + z₀ - ρ) / σ) ((0 - 1) / σ) ρ := fun ρ =>
    ((hasDerivAt_const ρ (S + z₀)).sub (hasDerivAt_id ρ)).div_const σ
  have hF : ∀ ρ, HasDerivAt F (F' ρ) ρ := fun ρ => by
    have h1 := ((hasDerivAt_pow 2 ρ).const_mul (n : ℝ)).add
      (((hasDerivAt_const ρ (S + z₀)).sub (hasDerivAt_id ρ)).pow 2)
    have h2 := ((hasDerivAt_robasBigOmega _).comp ρ (hin ρ)).const_mul T
    refine (h1.sub h2).congr_deriv ?_
    simp only [F', Pi.sub_apply, id, Nat.cast_ofNat, Nat.reduceSub, pow_one]
    rw [← hσ2]
    field_simp
    ring
  have hF' : ∀ ρ, HasDerivAt F' (2 * n + 2 - robasOmegaDeriv ((S + z₀ - ρ) / σ)) ρ := fun ρ => by
    have h1 := ((hasDerivAt_id ρ).const_mul (2 * (n : ℝ))).sub
      (((hasDerivAt_const ρ (S + z₀)).sub (hasDerivAt_id ρ)).const_mul 2)
    have h2 := ((hasDerivAt_robasOmega _).comp ρ (hin ρ)).const_mul σ
    refine (h1.add h2).congr_deriv ?_
    field_simp
    ring_nf
  have hmono : Monotone F' :=
    monotone_of_deriv_nonneg (fun ρ => (hF' ρ).differentiableAt) fun ρ => by
      rw [(hF' ρ).deriv]
      linarith [robasOmegaDeriv_le ((S + z₀ - ρ) / σ)]
  have hcrit : F' z₀ = 0 := by
    simp only [F', robasOmega]
    have e1 : (S + z₀ - z₀) / σ = S / σ := by ring
    have e2 : (S / σ) ^ 2 / 2 = S ^ 2 / (2 * T) := by rw [div_pow, hσ2]; field_simp
    rw [e1, e2]
    have : n * z₀ - S * robasG (S ^ 2 / (2 * T)) = 0 := hz₀
    field_simp
    linarith
  have key := le_of_monotone_deriv hF hmono hcrit (r j)
  have hΔ : robasDelta r j z₀ = F (r j) - F z₀ := by
    simp only [F, robasDelta]
    rw [← hT_def, ← hσ_def, ← hS]
    have e : (S + z₀ - z₀) / σ = S / σ := by ring
    rw [e]
    ring
  rw [hΔ]
  linarith

/-- **RoBAS-Full's conformal sets are intervals for `n ≥ 4`.** For calibration residuals not all
equal and every `α`, the set of candidate residuals `{z : p(z) > α}` is an interval, granted the
four certified bounds on the module's Dawson function (`RobasDawsonBounds`). Each unit's acceptance
set is an interval (`robas_acceptance_ordConnected`), they share a point (`robas_common_point`),
and Bersson and Hoff's Lemmas 1–2 finish. -/
theorem robasFull_interval (hb : RobasDawsonBounds) (hn : 4 ≤ n) (r : Fin n → ℝ)
    (hr : ∃ i k, r i ≠ r k) (α : ℝ) : Set.OrdConnected {z | α < robasPValue r z} := by
  obtain ⟨z₀, hz₀⟩ := robas_common_point r hr
  have hI := fun j => robas_acceptance_ordConnected hb hn r hr j
  have hnR : (n : ℝ) ≠ 0 := by
    have : (4 : ℝ) ≤ n := by exact_mod_cast hn
    positivity
  have hiff : ∀ z j, robasCandidateScore r z ≤ robasUnitScore r z j ↔ 0 ≤ robasDelta r j z := by
    intro z j
    have hT := robasSpread_pos r hr z
    rw [← robas_score_diff hnR r z hT j]
    constructor
    · intro h; exact mul_nonneg (by linarith) (sub_nonneg.2 h)
    · intro h
      have h2 : 0 < 2 * robasSpread r z := by linarith
      exact sub_nonneg.1 (le_of_mul_le_mul_left (by rw [mul_zero]; exact h) h2)
  have hm : (0 : ℝ) < n + 1 := by positivity
  have hcount : ∀ x w, (∀ j, 0 ≤ robasDelta r j x → 0 ≤ robasDelta r j w) →
      robasPValue r x ≤ robasPValue r w := by
    intro x w hxw
    unfold robasPValue
    refine div_le_div_of_nonneg_right ?_ hm.le
    have hsub : (Finset.univ.filter fun j => robasCandidateScore r x ≤ robasUnitScore r x j) ⊆
        (Finset.univ.filter fun j => robasCandidateScore r w ≤ robasUnitScore r w j) := by
      intro j hj
      simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hj ⊢
      rw [hiff] at hj ⊢
      exact hxw j hj
    have := Finset.card_le_card hsub
    have : ((Finset.univ.filter fun j => robasCandidateScore r x ≤ robasUnitScore r x j).card : ℝ) ≤
        (Finset.univ.filter fun j => robasCandidateScore r w ≤ robasUnitScore r w j).card := by
      exact_mod_cast this
    linarith
  refine ⟨fun x hx y hy w hw => ?_⟩
  simp only [Set.mem_setOf_eq] at hx hy ⊢
  rcases le_total w z₀ with h | h
  · exact lt_of_lt_of_le hx (hcount x w fun j hj => (hI j).out hj (hz₀ j) ⟨hw.1, h⟩)
  · exact lt_of_lt_of_le hy (hcount y w fun j hj => (hI j).out (hz₀ j) hj ⟨h, hw.2⟩)

end RobasFull

end Conformal
