import Enthymemata.TimeFrequency.ProfileBinaryLaws
import Enthymemata.TimeFrequency.HaarFamilies
import Mathlib.Combinatorics.SimpleGraph.Hasse

/-! # A path and its missing acoustic prediction -/

namespace TimeFrequency.LinkedProfiles
open MeasureTheory Set SimpleGraph
noncomputable section

def pathCenter (_ : Fin 3) : ℝ := 1 / 2
def flatContrast (_ : Fin 3) : ℝ := 1 / 16
def scaledContrast (v : Fin 3) : ℝ := if v = 1 then 1 / 32 else 1 / 8

private theorem flat_admissible (s : Bool) (v : Fin 3) :
    pathCenter v + sign s * flatContrast v ∈ Icc (1/4) (3/4) := by
  cases s <;> norm_num [Set.mem_Icc, pathCenter, flatContrast, sign]

private theorem scaled_admissible (s : Bool) (v : Fin 3) :
    pathCenter v + sign s * scaledContrast v ∈ Icc (1/4) (3/4) := by
  cases s <;> fin_cases v <;> norm_num [Set.mem_Icc, pathCenter, scaledContrast, sign]

theorem path_not_exchange : scaledContrast ≠ flatContrast ∧ scaledContrast ≠ -flatContrast := by
  constructor <;> intro h <;> have h₀ := congrFun h 0 <;>
    norm_num [flatContrast, scaledContrast] at h₀

theorem path_pair_laws {u v : Fin 3} (h : (pathGraph 3).Adj u v) :
    bitPairs (1/4) (3/4) pathCenter flatContrast u v =
      bitPairs (1/4) (3/4) pathCenter scaledContrast u v := by
  apply bitPairs_eq (by norm_num) _ _ _ flat_admissible scaled_admissible
  rw [pathGraph_adj] at h
  fin_cases u <;> fin_cases v <;> norm_num [flatContrast, scaledContrast] at *

def emittedBlock (bit : Bool) := Haar.signal .three (fraction (1/4) (3/4) bit) 0 0
def waveletReadout (bit : Bool) := Haar.normalized (emittedBlock bit)

theorem emitted_measurement (bit : Bool) :
    waveletReadout bit = Haar.features .three (fraction (1/4) (3/4) bit) := by
  apply Haar.normalized_signal
  cases bit <;> norm_num [Set.mem_Icc, fraction]

theorem waveletReadout_injective : Function.Injective waveletReadout := by
  intro i j he
  have ht (b : Bool) : fraction (1/4) (3/4) b ∈ Icc 0 1 := by
    cases b <;> norm_num [Set.mem_Icc, fraction]
  have hm : Haar.measured (emittedBlock i) = Haar.measured (emittedBlock j) := by
    simp only [emitted_measurement] at he
    simpa only [emittedBlock, Haar.measured_signal _ (ht _)] using he
  have h := (Haar.three_scale_identification _ _ (ht i) (ht j) 0 0 0 0 hm).2
  cases i <;> cases j <;> norm_num [fraction] at *

def acousticPair (a : Fin 3 → ℝ) (u v : Fin 3) :=
  readPairs waveletReadout (bitPairs (1/4) (3/4) pathCenter a u v)

theorem path_acoustic_laws {u v : Fin 3} (h : (pathGraph 3).Adj u v) :
    readPairs emittedBlock (bitPairs (1/4) (3/4) pathCenter flatContrast u v) =
      readPairs emittedBlock (bitPairs (1/4) (3/4) pathCenter scaledContrast u v) ∧
    acousticPair flatContrast u v = acousticPair scaledContrast u v := by
  simp only [acousticPair, path_pair_laws h, and_self]

theorem missing_edge_prediction :
    (acousticPair flatContrast 0 2).real {(waveletReadout true, waveletReadout true)} = 17/64 ∧
    (acousticPair scaledContrast 0 2).real {(waveletReadout true, waveletReadout true)} = 5/16 := by
  simp only [acousticPair, readPairs_cell _ waveletReadout_injective,
    bitPairs_cell (by norm_num : (1:ℝ)/4 < 3/4) _ _ flat_admissible,
    bitPairs_cell (by norm_num : (1:ℝ)/4 < 3/4) _ _ scaled_admissible]
  norm_num [success, bitMass, pathCenter, flatContrast, scaledContrast, sign,
    show (2 : Fin 3) ≠ 1 by decide]

/-- Closing 0–1–2 into a triangle identifies any admissible pair of source profiles. -/
theorem triangle_identification (P Q : ResponseModel (Fin 3)) (ha : ∀ v, P.contrast v ≠ 0)
    (same : ∀ u v, u ≠ v → Q.pairLaw u v = P.pairLaw u v) :
    Q.mean = P.mean ∨ Q.mean = fun s ↦ P.mean (!s) := by
  let G : SimpleGraph (Fin 3) := ⊤
  have e (i j : Fin 3) (hne : i ≠ j) : G.Adj i j := (top_adj ..).mpr hne
  let w : G.Walk 0 0 := .cons (e 0 1 (by decide))
    (.cons (e 1 2 (by decide)) (.cons (e 2 0 (by decide)) .nil))
  apply P.identify_profiles Q G preconnected_top (· ≠ ·) _ (fun _ _ h ↦ (top_adj ..).mp h)
    same ha 0 (Or.inr ⟨w, show Odd 3 from ⟨1, rfl⟩⟩)
  intro u
  fin_cases u
  · exact ⟨1, by decide⟩
  · exact ⟨2, by decide⟩
  · exact ⟨0, by decide⟩

/-- On the triangle, equality of the measured acoustic pair laws identifies the two
realized binary-response profiles, including when the candidate uses different centers. -/
theorem wavelet_triangle_profiles (μ a ν b : Fin 3 → ℝ)
    (hA : ∀ s v, μ v + sign s * a v ∈ Icc (1 / 4) (3 / 4))
    (hB : ∀ s v, ν v + sign s * b v ∈ Icc (1 / 4) (3 / 4)) (ha : ∀ v, a v ≠ 0)
    (same : ∀ i j : Fin 3, i ≠ j →
      readPairs waveletReadout (bitPairs (1/4) (3/4) ν b i j) =
        readPairs waveletReadout (bitPairs (1/4) (3/4) μ a i j)) :
    (∀ s v, ν v + sign s * b v = μ v + sign s * a v) ∨
      (∀ s v, ν v + sign s * b v = μ v + sign (!s) * a v) := by
  have hlh : (1:ℝ)/4 < 3/4 := by norm_num
  let P := binaryModel hlh μ a hA
  let Q := binaryModel hlh ν b hB
  have hcontrast (v : Fin 3) : P.contrast v ≠ 0 := by
    simpa [P, ResponseModel.contrast, binaryModel_mean, sign] using ha v
  have laws (i j : Fin 3) (hne : i ≠ j) : Q.pairLaw i j = P.pairLaw i j := by
    have hbits : bitPairs (1/4) (3/4) ν b i j = bitPairs (1/4) (3/4) μ a i j := by
      apply MeasureTheory.ext_iff_measureReal_singleton.mpr
      rintro ⟨x,y⟩
      rw [← readPairs_cell _ waveletReadout_injective _ x y, same i j hne,
        readPairs_cell _ waveletReadout_injective]
    simp only [P, Q, binaryModel_pairs, hbits]
  rcases triangle_identification P Q hcontrast laws with h | h
  · exact Or.inl fun s v ↦ by
      simpa only [P, Q, binaryModel_mean] using congrFun (congrFun h s) v
  · exact Or.inr fun s v ↦ by
      simpa only [P, Q, binaryModel_mean] using congrFun (congrFun h s) v

#print axioms wavelet_triangle_profiles
#print axioms path_not_exchange
#print axioms path_pair_laws
#print axioms emitted_measurement
#print axioms waveletReadout_injective
#print axioms path_acoustic_laws
#print axioms missing_edge_prediction
#print axioms triangle_identification

end
end TimeFrequency.LinkedProfiles
