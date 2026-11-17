import Enthymemata.TimeFrequency.MemoryDepthTree

/-! # Supported depths and simultaneous source-law confidence bands -/

namespace TimeFrequency.ContextTree.DepthEvidence

open Finset

open scoped Classical in
noncomputable def memoryDepth {D : ℕ} (K : Word D → Bool → ℝ) : ℕ :=
  (univ.filter (fun k : Fin (D + 1) ↦
    ∃ s ∈ minimalTree K, s.length = k.val)).sup Fin.val

theorem leaf_length_bound {D : ℕ} {K : Word D → Bool → ℝ} {s : List Bool}
    (hs : s ∈ minimalTree K) : s.length ≤ D := by
  obtain ⟨h, hh⟩ := hs.1
  simpa using hh.length_le

theorem leaf_length_le_depth {D : ℕ} {K : Word D → Bool → ℝ} {s : List Bool}
    (hs : s ∈ minimalTree K) : s.length ≤ memoryDepth K := by
  classical
  let k : Fin (D + 1) := ⟨s.length, by have hb := leaf_length_bound hs; omega⟩
  exact le_sup (f := Fin.val) (show k ∈ univ.filter (fun j : Fin (D + 1) ↦
    ∃ u ∈ minimalTree K, u.length = j.val) from mem_filter.mpr
      ⟨mem_univ _, s, hs, rfl⟩)

theorem memoryDepth_le_iff {D k : ℕ} (K : Word D → Bool → ℝ) :
    memoryDepth K ≤ k ↔ ∀ s ∈ minimalTree K, s.length ≤ k := by
  classical
  constructor
  · intro hk s hs
    exact (leaf_length_le_depth hs).trans hk
  · intro hk
    apply Finset.sup_le
    intro j hj
    obtain ⟨s, hs, he⟩ := (mem_filter.mp hj).2
    rw [← he]
    exact hk s hs

theorem memoryDepth_bound {D : ℕ} (K : Word D → Bool → ℝ) : memoryDepth K ≤ D :=
  (memoryDepth_le_iff K).mpr (fun _ hs ↦ leaf_length_bound hs)

def supportedDepths {Θ : Type*} (depth : Θ → ℕ) (C : Set Θ) : Set ℕ := depth '' C

private theorem supported_nonempty {Θ : Type*} (depth : Θ → ℕ) {C : Set Θ}
    (hC : C.Nonempty) : ∃ n, n ∈ supportedDepths depth C := by
  obtain ⟨θ, hθ⟩ := hC
  exact ⟨depth θ, θ, hθ, rfl⟩

open scoped Classical in
noncomputable def lowerSupported {Θ : Type*} (depth : Θ → ℕ) (C : Set Θ) : Option ℕ :=
  if hC : C.Nonempty then some (Nat.find (supported_nonempty depth hC)) else none

theorem lowerSupported_empty {Θ : Type*} (depth : Θ → ℕ) :
    lowerSupported depth (∅ : Set Θ) = none := by simp [lowerSupported]

theorem lowerSupported_spec {Θ : Type*} (depth : Θ → ℕ) (C : Set Θ) {k : ℕ}
    (hk : lowerSupported depth C = some k) :
    k ∈ supportedDepths depth C ∧ ∀ θ ∈ C, k ≤ depth θ := by
  classical
  unfold lowerSupported at hk
  split_ifs at hk with hC
  · have he := Option.some.inj hk
    rw [← he]
    exact ⟨Nat.find_spec (supported_nonempty depth hC),
      fun θ hθ ↦ Nat.find_min' _ ⟨θ, hθ, rfl⟩⟩

theorem covered_lower_depth {Θ : Type*} (depth : Θ → ℕ) {C : Set Θ} {θ : Θ}
    (hθ : θ ∈ C) : ∃ k, lowerSupported depth C = some k ∧ k ≤ depth θ := by
  classical
  have hC : C.Nonempty := ⟨θ, hθ⟩
  refine ⟨Nat.find (supported_nonempty depth hC), ?_, ?_⟩
  · simp [lowerSupported, hC]
  · exact Nat.find_min' _ ⟨θ, hθ, rfl⟩

def transitionBand {D : ℕ} (C : Set (Law (Word (D + 1)))) (qhat : Word D → ℝ)
    (t : ℝ) : Set (Law (Word (D + 1))) :=
  {P | P ∈ C ∧ ∀ h, |qhat h - conditional (historyResponseLaw P) h| < t}

def retainedSplit {D : ℕ} (C : Set (Law (Word (D + 1)))) (s : List Bool) : Prop :=
  C.Nonempty ∧ ∀ P ∈ C, ¬ constantSuffix (conditional (historyResponseLaw P)) s

theorem band_coverage {D : ℕ} (P : Law (Word (D + 1)))
    (C : Set (Law (Word (D + 1)))) (hP : P ∈ C) (n : ℕ) {t π δ : ℝ}
    (ht : 0 < t) (hu : t ≤ 1 / 2) (hπ : 0 < π)
    (hm : ∀ h, π ≤ historyMass (historyResponseLaw P) h) (hδ : 0 < δ)
    (hn : ⌈Real.log (2 * Fintype.card (Word D) / δ) / (π * t ^ 2)⌉₊ ≤ n) :
    probability (iid n (historyResponseLaw P))
      (fun x ↦ P ∉ transitionBand C (estimate x) t) ≤ δ := by
  have hb := transition_budget (historyResponseLaw P) n ht hu hπ hm hδ
    (Fintype.card_pos_iff.mpr inferInstance) hn
  apply (probability_mono _ ?_).trans hb
  intro x hx
  by_contra he
  apply hx
  exact ⟨hP, fun h ↦ lt_of_not_ge (fun hh ↦ he ⟨h, hh⟩)⟩

theorem band_lower_depth_error {D : ℕ} (P : Law (Word (D + 1)))
    (C : Set (Law (Word (D + 1)))) (n : ℕ) (t : ℝ) :
    probability (iid n (historyResponseLaw P)) (fun x ↦ ∃ k,
      lowerSupported (fun Q ↦ memoryDepth (nextKernel D Q.mass))
        (transitionBand C (estimate x) t) = some k ∧ memoryDepth (nextKernel D P.mass) < k) ≤
    probability (iid n (historyResponseLaw P))
      (fun x ↦ P ∉ transitionBand C (estimate x) t) := by
  apply probability_mono
  rintro x ⟨k, hk, hlarge⟩ hP
  exact (not_lt_of_ge ((lowerSupported_spec _ _ hk).2 P hP)) hlarge

theorem false_split_error {D : ℕ} (P : Law (Word (D + 1)))
    (C : Set (Law (Word (D + 1)))) (n : ℕ) (t : ℝ) (s : List Bool)
    (hc : constantSuffix (conditional (historyResponseLaw P)) s) :
    probability (iid n (historyResponseLaw P))
      (fun x ↦ retainedSplit (transitionBand C (estimate x) t) s) ≤
    probability (iid n (historyResponseLaw P))
      (fun x ↦ P ∉ transitionBand C (estimate x) t) := by
  apply probability_mono
  intro x hx hP
  exact hx.2 P hP hc

theorem band_split_power {D : ℕ} {P : Law (Word (D + 1))}
    {C : Set (Law (Word (D + 1)))} {qhat : Word D → ℝ} {t : ℝ} {s : List Bool}
    (hP : P ∈ transitionBand C qhat t)
    (hgap : ∃ h h' : Word D, s <:+ List.ofFn h ∧ s <:+ List.ofFn h' ∧
      4 * t < |conditional (historyResponseLaw P) h - conditional (historyResponseLaw P) h'|) :
    retainedSplit (transitionBand C qhat t) s := by
  refine ⟨⟨P, hP⟩, ?_⟩
  intro Q hQ hc
  obtain ⟨h, h', hh, hh', hg⟩ := hgap
  have he := hc h h' hh hh'
  have e1 := hP.2 h
  have e2 := hP.2 h'
  have e3 := hQ.2 h
  have e4 := hQ.2 h'
  have h1 := abs_sub_le (conditional (historyResponseLaw P) h) (qhat h)
    (conditional (historyResponseLaw P) h')
  have h2 := abs_sub_le (qhat h) (conditional (historyResponseLaw Q) h)
    (conditional (historyResponseLaw P) h')
  have h3 := abs_sub_le (conditional (historyResponseLaw Q) h') (qhat h')
    (conditional (historyResponseLaw P) h')
  rw [he] at h2 e3
  rw [abs_sub_comm (conditional (historyResponseLaw P) h) (qhat h)] at h1
  rw [abs_sub_comm (conditional (historyResponseLaw Q) h') (qhat h')] at h3
  linarith

theorem missed_split_error {D : ℕ} (P : Law (Word (D + 1)))
    (C : Set (Law (Word (D + 1)))) (n : ℕ) (t : ℝ) (s : List Bool)
    (hg : ∃ h h' : Word D, s <:+ List.ofFn h ∧ s <:+ List.ofFn h' ∧
      4 * t < |conditional (historyResponseLaw P) h - conditional (historyResponseLaw P) h'|) :
    probability (iid n (historyResponseLaw P))
      (fun x ↦ ¬ retainedSplit (transitionBand C (estimate x) t) s) ≤
    probability (iid n (historyResponseLaw P))
      (fun x ↦ P ∉ transitionBand C (estimate x) t) := by
  apply probability_mono
  intro x hx hP
  exact hx (band_split_power hP hg)

theorem common_law_depth_ambiguity {X : Type*} [Fintype X] (R : Law X)
    (S : X → Set ℕ) (d : ℕ) {δ : ℝ}
    (h0 : probability R (fun x ↦ 0 ∉ S x) ≤ δ)
    (hd : probability R (fun x ↦ d ∉ S x) ≤ δ) :
    probability R (fun x ↦ ¬(0 ∈ S x ∧ d ∈ S x)) ≤ 2 * δ := by
  classical
  have he : (fun x ↦ ¬(0 ∈ S x ∧ d ∈ S x)) =
      (fun x ↦ 0 ∉ S x ∨ d ∉ S x) := by
    funext x
    exact propext (not_and_or)
  rw [he]
  exact (probability_union R _ _).trans (by linarith)

#print axioms missed_split_error
#print axioms common_law_depth_ambiguity

#print axioms leaf_length_bound
#print axioms leaf_length_le_depth
#print axioms memoryDepth_le_iff
#print axioms memoryDepth_bound
#print axioms lowerSupported_empty
#print axioms lowerSupported_spec
#print axioms covered_lower_depth
#print axioms band_coverage
#print axioms band_lower_depth_error
#print axioms false_split_error
#print axioms band_split_power

end TimeFrequency.ContextTree.DepthEvidence
