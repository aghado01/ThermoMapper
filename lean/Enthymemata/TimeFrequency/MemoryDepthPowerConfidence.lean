import Enthymemata.TimeFrequency.MemoryDepthStability
import Enthymemata.TimeFrequency.MemoryDepthUnpaired
import Enthymemata.TimeFrequency.MemoryDepthPrediction

/-! # Valid source confidence sets from the unpaired depth-four experiment -/

namespace TimeFrequency.ContextTree.DepthEvidence

def powerBand (C : Set Source) (m : Fin 8 → ℝ) (ε : ℝ) : Set Source :=
  {P | P ∈ C ∧ ∀ j, |m j - powerSummaries P j| < ε}

def retainedPowerSplit (C : Set Source) (s : List Bool) : Prop :=
  C.Nonempty ∧ ∀ P ∈ C, ¬ constantSuffix (fun h ↦ nextKernel 4 (P.law 5) h true) s

theorem power_band_coverage (P : Source) (hp : ∀ x, 0 < P.law 5 x)
    (C : Set Source) (hP : P ∈ C) (n : ℕ) {ε δ : ℝ}
    (he : 0 < ε) (hδ : 0 < δ) (hδu : δ < 1)
    (hn : ⌈Real.log (16 / δ) / (2 * ε ^ 2)⌉₊ ≤ n) :
    probability (iid 6 (iid n (sourceBlockLaw P hp)))
      (fun x ↦ P ∉ powerBand C (recordedSummaries x) ε) ≤ δ := by
  apply (probability_mono _ ?_).trans (unpaired_summary_budget P hp n he hδ hδu hn)
  intro x hx
  by_contra hb
  exact hx ⟨hP, fun j ↦ lt_of_not_ge (fun hj ↦ hb ⟨j, hj⟩)⟩

theorem power_band_kernel_diameter {P Q : Source} {C : Set Source}
    {m : Fin 8 → ℝ} {η ε : ℝ} (hη : 0 < η) (he : 0 ≤ ε)
    (hP : P ∈ powerBand C m ε) (hQ : Q ∈ powerBand C m ε)
    (hp : ∀ h : Word 4, η ≤ P.law 4 h) (hq : ∀ h : Word 4, η ≤ Q.law 4 h)
    (h : Word 4) :
    |nextKernel 4 (P.law 5) h true - nextKernel 4 (Q.law 5) h true| ≤
      2 * (385 / 16 : ℝ) * ε / η := by
  have hm (j : Fin 8) : |powerSummaries P j - powerSummaries Q j| ≤ 2 * ε := by
    have ht := abs_sub_le (powerSummaries P j) (m j) (powerSummaries Q j)
    rw [abs_sub_comm (powerSummaries P j) (m j)] at ht
    linarith [hP.2 j, hQ.2 j]
  have hb := source_kernel_stability P Q hη (show 0 ≤ 2 * ε by positivity) hp hq hm h
  convert hb using 1
  ring

theorem power_band_lower_depth_error (P : Source) (hp : ∀ x, 0 < P.law 5 x)
    (C : Set Source) (n : ℕ) (ε : ℝ) :
    probability (iid 6 (iid n (sourceBlockLaw P hp))) (fun x ↦ ∃ k,
      lowerSupported (fun Q ↦ memoryDepth (nextKernel 4 (Q.law 5)))
        (powerBand C (recordedSummaries x) ε) = some k ∧
          memoryDepth (nextKernel 4 (P.law 5)) < k) ≤
    probability (iid 6 (iid n (sourceBlockLaw P hp)))
      (fun x ↦ P ∉ powerBand C (recordedSummaries x) ε) := by
  apply probability_mono
  rintro x ⟨k, hk, hlarge⟩ hP
  exact (not_lt_of_ge ((lowerSupported_spec _ _ hk).2 P hP)) hlarge

theorem power_false_split_error (P : Source) (hp : ∀ x, 0 < P.law 5 x)
    (C : Set Source) (n : ℕ) (ε : ℝ) (s : List Bool)
    (hc : constantSuffix (fun h ↦ nextKernel 4 (P.law 5) h true) s) :
    probability (iid 6 (iid n (sourceBlockLaw P hp)))
      (fun x ↦ retainedPowerSplit (powerBand C (recordedSummaries x) ε) s) ≤
    probability (iid 6 (iid n (sourceBlockLaw P hp)))
      (fun x ↦ P ∉ powerBand C (recordedSummaries x) ε) := by
  apply probability_mono
  intro x hx hP
  exact hx.2 P hP hc

theorem power_split_power {P : Source} {C : Set Source} {m : Fin 8 → ℝ} {η ε : ℝ}
    (hη : 0 < η) (he : 0 ≤ ε) (hf : ∀ Q ∈ C, ∀ h : Word 4, η ≤ Q.law 4 h)
    {s : List Bool} (hP : P ∈ powerBand C m ε)
    (hg : ∃ h h' : Word 4, s <:+ List.ofFn h ∧ s <:+ List.ofFn h' ∧
      4 * (385 / 16 : ℝ) * ε / η <
        |nextKernel 4 (P.law 5) h true - nextKernel 4 (P.law 5) h' true|) :
    retainedPowerSplit (powerBand C m ε) s := by
  refine ⟨⟨P, hP⟩, ?_⟩
  intro Q hQ hc
  obtain ⟨h, h', hh, hh', hg⟩ := hg
  have h1 := power_band_kernel_diameter hη he hP hQ (hf P hP.1) (hf Q hQ.1) h
  have h2 := power_band_kernel_diameter hη he hP hQ (hf P hP.1) (hf Q hQ.1) h'
  have ht := abs_sub_le (nextKernel 4 (P.law 5) h true)
    (nextKernel 4 (Q.law 5) h true) (nextKernel 4 (P.law 5) h' true)
  have heq : nextKernel 4 (Q.law 5) h true = nextKernel 4 (Q.law 5) h' true :=
    hc h h' hh hh'
  rw [heq] at h1
  rw [heq, abs_sub_comm (nextKernel 4 (Q.law 5) h' true)
    (nextKernel 4 (P.law 5) h' true)] at ht
  have hb : |nextKernel 4 (P.law 5) h true - nextKernel 4 (P.law 5) h' true| ≤
      4 * (385 / 16 : ℝ) * ε / η := calc
    _ ≤ _ := ht
    _ ≤ 2 * (385 / 16 : ℝ) * ε / η + 2 * (385 / 16 : ℝ) * ε / η :=
      add_le_add h1 h2
    _ = _ := by ring
  exact (not_lt_of_ge hb) hg

theorem power_missed_split_error (P : Source) (hp : ∀ x, 0 < P.law 5 x)
    (C : Set Source) (n : ℕ) {η ε : ℝ} (hη : 0 < η) (he : 0 ≤ ε)
    (hf : ∀ Q ∈ C, ∀ h : Word 4, η ≤ Q.law 4 h) (s : List Bool)
    (hg : ∃ h h' : Word 4, s <:+ List.ofFn h ∧ s <:+ List.ofFn h' ∧
      4 * (385 / 16 : ℝ) * ε / η <
        |nextKernel 4 (P.law 5) h true - nextKernel 4 (P.law 5) h' true|) :
    probability (iid 6 (iid n (sourceBlockLaw P hp)))
      (fun x ↦ ¬ retainedPowerSplit (powerBand C (recordedSummaries x) ε) s) ≤
    probability (iid 6 (iid n (sourceBlockLaw P hp)))
      (fun x ↦ P ∉ powerBand C (recordedSummaries x) ε) := by
  apply probability_mono
  intro x hx hP
  exact hx (power_split_power hη he hf hP hg)

theorem source_binary_kernel (P : Source) {η : ℝ} (hη : 0 < η)
    (hp : ∀ h : Word 4, η ≤ P.law 4 h) :
    binaryKernel (fun h ↦ nextKernel 4 (P.law 5) h true) = nextKernel 4 (P.law 5) := by
  funext h b
  have hd : 0 < P.law 5 (Fin.snoc h false) + P.law 5 (Fin.snoc h true) := by
    rw [prefix_mass]
    exact hη.trans_le (hp h)
  cases b
  · change 1 - P.law 5 (Fin.snoc h true) /
        (P.law 5 (Fin.snoc h false) + P.law 5 (Fin.snoc h true)) =
      P.law 5 (Fin.snoc h false) /
        (P.law 5 (Fin.snoc h false) + P.law 5 (Fin.snoc h true))
    field_simp [hd.ne']
    ring
  · rfl

theorem power_candidate_tree {P Q : Source} {C : Set Source} {m : Fin 8 → ℝ}
    {η ε γ : ℝ} (hη : 0 < η) (he : 0 ≤ ε)
    (hP : P ∈ powerBand C m ε) (hQ : Q ∈ powerBand C m ε)
    (hp : ∀ h : Word 4, η ≤ P.law 4 h) (hq : ∀ h : Word 4, η ≤ Q.law 4 h)
    (hs : separated (fun h ↦ nextKernel 4 (P.law 5) h true) γ)
    (hm : 2 * (385 / 16 : ℝ) * ε / η < γ / 4) :
    suffixTree (d := 4)
      (empiricalConstant (fun h ↦ nextKernel 4 (Q.law 5) h true) γ) =
        minimalTree (nextKernel 4 (P.law 5)) := by
  have herr (h : Word 4) :
      |nextKernel 4 (Q.law 5) h true - nextKernel 4 (P.law 5) h true| < γ / 4 := by
    rw [abs_sub_comm]
    exact (power_band_kernel_diameter hη he hP hQ hp hq h).trans_lt hm
  have ht := threshold_tree hs herr
  rwa [source_binary_kernel P hη hp] at ht

theorem power_candidate_prediction {P Q : Source} {C : Set Source} {m : Fin 8 → ℝ}
    {η ε : ℝ} (hη : 0 < η) (he : 0 ≤ ε)
    (hP : P ∈ powerBand C m ε) (hQ : Q ∈ powerBand C m ε)
    (hp : ∀ h : Word 4, η ≤ P.law 4 h) (hq : ∀ h : Word 4, η ≤ Q.law 4 h) :
    sourceExcess P (fun h ↦ nextKernel 4 (Q.law 5) h true) ≤
      (2 * (385 / 16 : ℝ) * ε / η) ^ 2 := by
  apply uniform_prediction_tolerance
  intro h
  rw [abs_sub_comm]
  exact power_band_kernel_diameter hη he hP hQ hp hq h

#print axioms source_binary_kernel
#print axioms power_candidate_tree
#print axioms power_candidate_prediction

#print axioms power_band_coverage
#print axioms power_band_kernel_diameter
#print axioms power_band_lower_depth_error
#print axioms power_false_split_error
#print axioms power_split_power
#print axioms power_missed_split_error

end TimeFrequency.ContextTree.DepthEvidence
