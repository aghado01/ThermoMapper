import Enthymemata.TimeFrequency.MemoryDepthUpper
import Enthymemata.TimeFrequency.MemoryDepthObservations

/-! # Observation budgets for detecting the stationary depth-d alternative -/

namespace TimeFrequency.ContextTree.DepthEvidence

open Finset Collision

variable {X : Type*} [Fintype X]

theorem common_only_obstruction (R : Law X) (t : X → ℝ) {δ : ℝ} (hδ : δ < 1 / 2) :
    ¬(expect R t ≤ δ ∧ expect R (fun x ↦ 1 - t x) ≤ δ) := by
  have he : expect R t + expect R (fun x ↦ 1 - t x) = 1 := by
    simp only [expect, mul_sub, mul_one, sum_sub_distrib, R.total]
    ring
  rintro ⟨h0, h1⟩
  linarith

open scoped Classical in
/-- A fixed separate-view allocation. Each rectangular vector and each taper bit comes
from a fresh stationary block. The rectangular records form a common ancillary law. -/
theorem physical_budgets {m : ℕ} {a δ : ℝ} (ha : 0 < a) (hu : a < 1 / 2)
    (hδ : 0 < δ) (hδu : δ < 1 / 2) (r n : ℕ) (w o : Fin r → ℕ)
    (ho : ∀ i, o i + w i ≤ L m) :
    let R := independent (fun i ↦ rectangularLaw (fairLaw (L m)) (w i) (o i) (ho i))
    (∀ t, (∀ x, t x ∈ Set.Icc (0 : ℝ) 1) →
      1 - δ ≤ expect
        (product R (iid n (bernoulli ((1 + a) / 8) (taper_parameters ha hu).2))) t →
      expect (product R (iid n (bernoulli (1 / 8) (taper_parameters ha hu).1))) t ≤ δ →
      ⌈7 * kl (1 - δ) δ / a ^ 2⌉₊ ≤ n) ∧
    (⌈128 * Real.log (1 / δ) / a ^ 2⌉₊ ≤ n →
      expect (product R (iid n (bernoulli (1 / 8) (taper_parameters ha hu).1)))
          (fun z ↦ taperDecision n a z.2) ≤ δ ∧
      expect (product R (iid n (bernoulli ((1 + a) / 8) (taper_parameters ha hu).2)))
          (fun z ↦ 1 - taperDecision n a z.2) ≤ δ) := by
  dsimp only
  constructor
  · intro t ht hP hQ
    exact necessary_budget_rounded ha hu hδ hδu n _ t ht hP hQ
  · intro hn
    exact sufficient_budget_rounded ha hu hδ n _ hn

#print axioms common_only_obstruction
#print axioms physical_budgets

end TimeFrequency.ContextTree.DepthEvidence
