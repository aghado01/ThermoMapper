import Mathlib.Analysis.Convex.PathConnected
import Mathlib.Analysis.Normed.Module.FiniteDimension
import Mathlib.Combinatorics.SetFamily.Shatter
import Mathlib.LinearAlgebra.Dual.Lemmas
import Mathlib.Topology.Order.IntermediateValue

/-! # Strict signs of a central hyperplane arrangement (S22) -/

open Finset Module

namespace GenerativePrior

variable {ι E : Type*} [Fintype ι] [DecidableEq ι]
  [NormedAddCommGroup E] [NormedSpace ℝ E]
  (A : ι → E →ₗ[ℝ] ℝ)

/-- Indices in the mask have strictly positive response, all others strictly negative. -/
def signCell (s : Finset ι) : Set E :=
  {x | ∀ i, if i ∈ s then 0 < A i x else A i x < 0}

noncomputable def strictPatterns : Finset (Finset ι) := by
  classical exact univ.filter fun s ↦ (signCell A s).Nonempty

theorem mem_strictPatterns {s : Finset ι} :
    s ∈ strictPatterns A ↔ (signCell A s).Nonempty := by classical simp [strictPatterns]

theorem compl_mem_strictPatterns {s : Finset ι} (hs : s ∈ strictPatterns A) :
    sᶜ ∈ strictPatterns A := by
  obtain ⟨x, hx⟩ := (mem_strictPatterns A).mp hs
  refine (mem_strictPatterns A).mpr ⟨-x, fun i ↦ ?_⟩
  specialize hx i
  by_cases hi : i ∈ s <;> simpa [hi] using hx

/-- A shattered set of row indices supports no nonzero linear relation. -/
theorem independent_of_shatters {s : Finset ι} (hs : (strictPatterns A).Shatters s) :
    LinearIndependent ℝ (fun i : s ↦ A i.val) := by
  classical
  apply Fintype.linearIndependent_iff.mpr
  intro c hc
  let coeff (j : ι) := if h : j ∈ s then c ⟨j, h⟩ else 0
  obtain ⟨u, hu, htrace⟩ := hs (filter_subset (fun j ↦ 0 < coeff j) s)
  obtain ⟨x, hx⟩ := (mem_strictPatterns A).mp hu
  have hsign (i : s) : i.val ∈ u ↔ 0 < c i := by
    have h := Finset.ext_iff.mp htrace i.val
    simpa [i.property, coeff] using h
  have hzero : ∑ i : s, c i * A i.val x = 0 := by
    simpa [LinearMap.sum_apply] using congrArg (fun f : E →ₗ[ℝ] ℝ ↦ f x) hc
  have hterms (i : s) : 0 ≤ c i * A i.val x := by
    have hi := hx i.val
    by_cases hp : 0 < c i
    · exact (mul_pos hp (by simpa [(hsign i).mpr hp] using hi)).le
    · exact mul_nonneg_of_nonpos_of_nonpos (le_of_not_gt hp)
        (le_of_lt (by simpa [mt (hsign i).mp hp] using hi))
  intro i
  have hz := (sum_eq_zero_iff_of_nonneg (fun j _ ↦ hterms j)).mp hzero i (mem_univ i)
  have hi := hx i.val
  apply (mul_eq_zero.mp hz).resolve_right
  by_cases hp : i.val ∈ u
  · exact ne_of_gt (by simpa [hp] using hi)
  · exact ne_of_lt (by simpa [hp] using hi)

/-- Positive masks for a distinguished first hyperplane. -/
noncomputable def positivePatterns (a : ι) := (strictPatterns A).filter (a ∈ ·)

theorem card_strictPatterns (a : ι) :
    (strictPatterns A).card = 2 * (positivePatterns A a).card := by
  classical
  have hpair : ((strictPatterns A).filter (a ∉ ·)).card = (positivePatterns A a).card := by
    apply card_bij (fun s _ ↦ sᶜ)
    · intro s hs
      obtain ⟨hs, ha⟩ := mem_filter.mp hs
      exact mem_filter.mpr ⟨compl_mem_strictPatterns A hs, mem_compl.mpr ha⟩
    · intro s _ t _ he
      simpa using congrArg (fun u : Finset ι ↦ uᶜ) he
    · intro s hs
      obtain ⟨hs, ha⟩ := mem_filter.mp hs
      exact ⟨sᶜ, mem_filter.mpr ⟨compl_mem_strictPatterns A hs, by simpa using ha⟩, compl_compl s⟩
  have h := card_filter_add_card_filter_not (s := strictPatterns A) (p := (a ∈ ·))
  change (positivePatterns A a).card + _ = _ at h
  omega

/-- Antipodal symmetry adds the distinguished row to every shattered trace of the positive masks. -/
theorem shatters_insert_first {a : ι} {s : Finset ι} (hs : (positivePatterns A a).Shatters s) :
    a ∉ s ∧ (strictPatterns A).Shatters (insert a s) := by
  classical
  obtain ⟨u, hu, he⟩ := hs (empty_subset s)
  have ha := (mem_filter.mp hu).2
  have hnot : a ∉ s := fun h ↦ by
    have hboth := mem_inter.mpr ⟨h, ha⟩
    simp [he] at hboth
  refine ⟨hnot, ?_⟩
  intro t ht
  by_cases hat : a ∈ t
  · obtain ⟨v, hv, he⟩ := hs (show s ∩ t ⊆ s from inter_subset_left)
    obtain ⟨hv, hav⟩ := mem_filter.mp hv
    refine ⟨v, hv, ?_⟩
    ext i
    have hi := Finset.ext_iff.mp he i
    have hit := @ht i
    simp only [mem_inter, mem_insert] at *
    grind
  · obtain ⟨v, hv, he⟩ := hs (Finset.sdiff_subset (s := s) (t := t))
    obtain ⟨hv, hav⟩ := mem_filter.mp hv
    refine ⟨vᶜ, compl_mem_strictPatterns A hv, ?_⟩
    ext i
    have hi := Finset.ext_iff.mp he i
    have hit := @ht i
    simp only [mem_inter, mem_sdiff, mem_compl, mem_insert] at *
    grind

theorem strictPatterns_bound [FiniteDimensional ℝ E] (a : ι) :
    (strictPatterns A).card ≤ 2 * ∑ j ∈ range (finrank ℝ E), (Fintype.card ι - 1).choose j := by
  classical
  rw [card_strictPatterns A a]
  apply Nat.mul_le_mul_left
  let tail := (univ : Finset ι).erase a
  have htail : tail.card = Fintype.card ι - 1 := by simp [tail]
  simp_rw [← htail, ← card_powersetCard]
  refine ((positivePatterns A a).card_le_card_shatterer.trans (card_le_card ?_)).trans
    (card_biUnion_le (s := range (finrank ℝ E)) (t := fun j ↦ tail.powersetCard j))
  intro s hs
  obtain ⟨ha, hfull⟩ := shatters_insert_first A (mem_shatterer.mp hs)
  have hdim := (independent_of_shatters A hfull).fintype_card_le_finrank
  change Fintype.card (insert a s : Finset ι) ≤ finrank ℝ (Module.Dual ℝ E) at hdim
  rw [Subspace.dual_finrank_eq] at hdim
  rw [Fintype.card_coe, card_insert_of_notMem ha] at hdim
  refine mem_biUnion.mpr ⟨s.card, mem_range.mpr (by omega), mem_powersetCard.mpr ⟨?_, rfl⟩⟩
  intro i hi
  exact mem_erase.mpr ⟨fun he ↦ ha (he ▸ hi), mem_univ _⟩

/-- The complement of the union of the kernels. -/
def offHyperplanes : Set E := {x | ∀ i, A i x ≠ 0}

omit [Fintype ι] in
theorem signCell_subset_offHyperplanes (s : Finset ι) : signCell A s ⊆ offHyperplanes A := by
  intro x hx i
  have hi := hx i
  by_cases h : i ∈ s
  · exact ne_of_gt (by simpa [h] using hi)
  · exact ne_of_lt (by simpa [h] using hi)

omit [Fintype ι] in
theorem convex_signCell (s : Finset ι) : Convex ℝ (signCell A s) := by
  unfold signCell
  rw [Set.setOf_forall]
  apply convex_iInter
  intro i
  by_cases h : i ∈ s
  · simpa only [if_pos h] using convex_halfSpace_gt (A i).isLinear 0
  · simpa only [if_neg h] using convex_halfSpace_lt (A i).isLinear 0

section Chambers

variable [FiniteDimensional ℝ E]

omit [Fintype ι] in
theorem isOpen_signCell [Finite ι] (s : Finset ι) : IsOpen (signCell A s) := by
  unfold signCell
  rw [Set.setOf_forall]
  apply isOpen_iInter_of_finite
  intro i
  have hc := (A i).continuous_of_finiteDimensional
  by_cases h : i ∈ s
  · simpa only [if_pos h] using isOpen_lt continuous_const hc
  · simpa only [if_neg h] using isOpen_lt hc continuous_const

omit [Fintype ι] in
theorem component_eq_signCell {s : Finset ι} {x : E} (hx : x ∈ signCell A s) :
    connectedComponentIn (offHyperplanes A) x = signCell A s := by
  let C := connectedComponentIn (offHyperplanes A) x
  have hxC : x ∈ C := mem_connectedComponentIn (signCell_subset_offHyperplanes A s hx)
  have hsub : C ⊆ offHyperplanes A := connectedComponentIn_subset _ _
  have hconn : IsPreconnected C := isPreconnected_connectedComponentIn
  refine Set.Subset.antisymm ?_ ((convex_signCell A s).isPreconnected.subset_connectedComponentIn
    hx (signCell_subset_offHyperplanes A s))
  intro y hy i
  have hcross {p q : E} (hp : p ∈ C) (hq : q ∈ C) (hpq : 0 ∈ Set.Icc (A i p) (A i q)) :
      False := by
    obtain ⟨z, hz, hzero⟩ := hconn.intermediate_value hp hq
      (A i).continuous_of_finiteDimensional.continuousOn hpq
    exact hsub hz i hzero
  have hi := hx i
  by_cases ha : i ∈ s
  · rw [if_pos ha] at hi ⊢
    exact lt_of_not_ge fun hn ↦ hcross hy hxC ⟨hn, hi.le⟩
  · rw [if_neg ha] at hi ⊢
    exact lt_of_not_ge fun hn ↦ hcross hxC hy ⟨hi.le, hn⟩

/-- Enumerate the nonempty cells as subsets of the ambient space. -/
noncomputable def signRegions : Finset (Set E) := by
  classical exact (strictPatterns A).image (signCell A)

theorem mem_signRegions {R : Set E} :
    R ∈ signRegions A ↔ ∃ x ∈ offHyperplanes A, connectedComponentIn (offHyperplanes A) x = R := by
  classical
  constructor
  · intro hR
    obtain ⟨s, hs, rfl⟩ := mem_image.mp hR
    obtain ⟨x, hx⟩ := (mem_strictPatterns A).mp hs
    exact ⟨x, signCell_subset_offHyperplanes A s hx, component_eq_signCell A hx⟩
  · rintro ⟨x, hx, rfl⟩
    let s := univ.filter fun i ↦ 0 < A i x
    have hcell : x ∈ signCell A s := by
      intro i
      by_cases hi : 0 < A i x
      · simp [s, hi]
      · simpa [s, hi] using lt_of_le_of_ne (le_of_not_gt hi) (hx i)
    exact mem_image.mpr ⟨s, (mem_strictPatterns A).mpr ⟨x, hcell⟩,
      (component_eq_signCell A hcell).symm⟩

omit [FiniteDimensional ℝ E] in
theorem card_signRegions : (signRegions A).card = (strictPatterns A).card := by
  classical
  apply card_image_of_injOn
  intro s hs t _ he
  obtain ⟨x, hx⟩ := (mem_strictPatterns A).mp hs
  have hy : x ∈ signCell A t := he ▸ hx
  ext i
  have h := hx i
  have h' := hy i
  by_cases hi : i ∈ s <;> by_cases hj : i ∈ t <;> simp_all
    <;> linarith

/-- The central-arrangement bound with arbitrary repetitions and linear dependencies. -/
theorem region_bound (a : ι) :
    (signRegions A).card ≤ 2 * ∑ j ∈ range (finrank ℝ E), (Fintype.card ι - 1).choose j :=
  (card_signRegions A).le.trans (strictPatterns_bound A a)

omit [FiniteDimensional ℝ E] in
theorem empty_arrangement [IsEmpty ι] : signRegions A = {Set.univ} := by
  classical
  have hp : strictPatterns A = {∅} := by
    ext s
    have hs : s = ∅ := Subsingleton.elim _ _
    simp only [hs, mem_strictPatterns, mem_singleton, iff_true]
    exact ⟨0, isEmptyElim⟩
  have hc : signCell A ∅ = Set.univ := by ext; simp [signCell]
  simp [signRegions, hp, hc]

theorem zero_dimensional_arrangement (a : ι) (hk : finrank ℝ E = 0) :
    signRegions A = ∅ := by
  have h := region_bound A a
  simpa [hk] using h

theorem no_regions_of_zero_row {a : ι} (ha : A a = 0) : signRegions A = ∅ := by
  apply eq_empty_iff_forall_notMem.mpr
  intro R hR
  obtain ⟨x, hx, _⟩ := (mem_signRegions A).mp hR
  exact hx a (by simp [ha])

end Chambers

omit [Fintype ι] in
/-- On a strict cell the bias-free ReLU layer agrees with a single linear map. -/
theorem relu_eq_masked_linear {s : Finset ι} {x : E} (hx : x ∈ signCell A s) :
    (fun i ↦ max 0 (A i x)) = (LinearMap.pi fun i ↦ if i ∈ s then A i else 0) x := by
  funext i
  have h := hx i
  change max 0 (A i x) = (if i ∈ s then A i else 0) x
  split_ifs with hi
  · rw [if_pos hi] at h
    exact max_eq_right h.le
  · rw [if_neg hi] at h
    exact max_eq_left h.le

/-- Coordinate hyperplanes realize all orthants, including the four quadrants in the plane. -/
theorem coordinate_regions :
    (signRegions (fun i : ι ↦ LinearMap.proj i : ι → (ι → ℝ) →ₗ[ℝ] ℝ)).card =
      2 ^ Fintype.card ι := by
  rw [card_signRegions]
  have h : strictPatterns (fun i : ι ↦ LinearMap.proj i : ι → (ι → ℝ) →ₗ[ℝ] ℝ) = univ := by
    apply eq_univ_of_forall
    intro s
    refine (mem_strictPatterns _).mpr ⟨fun i ↦ if i ∈ s then 1 else -1, fun i ↦ ?_⟩
    by_cases hi : i ∈ s <;> simp [hi]
  simp [h]

#print axioms mem_strictPatterns
#print axioms compl_mem_strictPatterns
#print axioms independent_of_shatters
#print axioms card_strictPatterns
#print axioms shatters_insert_first
#print axioms strictPatterns_bound
#print axioms signCell_subset_offHyperplanes
#print axioms convex_signCell
#print axioms isOpen_signCell
#print axioms component_eq_signCell
#print axioms mem_signRegions
#print axioms card_signRegions
#print axioms region_bound
#print axioms empty_arrangement
#print axioms zero_dimensional_arrangement
#print axioms no_regions_of_zero_row
#print axioms relu_eq_masked_linear
#print axioms coordinate_regions

end GenerativePrior
