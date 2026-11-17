import Enthymemata.TimeFrequency.MemoryDepthLikelihood

/-! # Finite joint activity laws and the three observable targets

The future alphabet is ordered 00, 10, 01, 11. Zero future probabilities
are allowed. Strict positivity applies only to the retained history support.
-/

namespace TimeFrequency.PredictiveActivity

open Finset hiding expect
open ContextTree.DepthEvidence

structure ActivityLaw where
  prob : Fin 4 → ℝ
  nonneg : ∀ y, 0 ≤ prob y
  total : ∑ y, prob y = 1

noncomputable def targetMeans (p : ActivityLaw) : Fin 3 → ℝ :=
  ![p.prob 1 + p.prob 3, p.prob 2 + p.prob 3, p.prob 3]

def invertTargets {R : Type*} [Ring R] (m : Fin 3 → R) : Fin 4 → R :=
  ![1 - m 0 - m 1 + m 2, m 0 - m 2, m 1 - m 2, m 2]

theorem targets_determine_law (p : ActivityLaw) :
    invertTargets (targetMeans p) = p.prob := by
  have ht := p.total
  norm_num [Fin.sum_univ_succ] at ht
  change p.prob 0 + (p.prob 1 + (p.prob 2 + p.prob 3)) = 1 at ht
  funext y
  fin_cases y <;> simp [invertTargets, targetMeans]
  linarith

theorem three_target_injective (p q : ActivityLaw)
    (h : targetMeans p = targetMeans q) : p = q := by
  have he : p.prob = q.prob := by rw [← targets_determine_law p, h, targets_determine_law]
  cases p
  cases q
  cases he
  rfl

theorem inversion_squared_bound (a b : Fin 3 → ℝ) :
    (∑ y, (invertTargets a y - invertTargets b y) ^ 2) ≤
      8 * ∑ j, (a j - b j) ^ 2 := by
  norm_num [Fin.sum_univ_succ, invertTargets]
  nlinarith [sq_nonneg (a 0 - b 0), sq_nonneg (a 1 - b 1),
    sq_nonneg ((a 0 - b 0) - (a 1 - b 1)),
    sq_nonneg ((a 0 - b 0) + (a 2 - b 2)),
    sq_nonneg ((a 1 - b 1) + (a 2 - b 2))]

theorem activity_squared_bound (p q : ActivityLaw) :
    (∑ y, (p.prob y - q.prob y) ^ 2) ≤
      8 * ∑ j, (targetMeans p j - targetMeans q j) ^ 2 := by
  have h := inversion_squared_bound (targetMeans p) (targetMeans q)
  rwa [targets_determine_law, targets_determine_law] at h

def binaryRisk (u a : ℝ) : ℝ := u * (1 - a) ^ 2 + (1 - u) * a ^ 2

theorem binary_excess (u a : ℝ) :
    binaryRisk u a - binaryRisk u u = (a - u) ^ 2 := by
  unfold binaryRisk
  ring

def outcomeA (y : Fin 4) : Bool := decide (y = 1 ∨ y = 3)
def outcomeB (y : Fin 4) : Bool := decide (y = 2 ∨ y = 3)

def observableTarget (j : Fin 3) (y : Fin 4) : ℝ :=
  if (if j = 0 then outcomeA y else if j = 1 then outcomeB y
    else outcomeA y && outcomeB y) then 1 else 0

private theorem sum_four (f : Fin 4 → ℝ) :
    (∑ y, f y) = f 0 + f 1 + f 2 + f 3 := by
  simp only [Fin.sum_univ_succ, Fin.sum_univ_zero, add_zero]
  change f 0 + (f 1 + (f 2 + f 3)) = _
  ring

theorem target_means_observable (p : ActivityLaw) (j : Fin 3) :
    (∑ y, p.prob y * observableTarget j y) = targetMeans p j := by
  rw [sum_four]
  fin_cases j <;> norm_num [observableTarget, outcomeA, outcomeB, targetMeans, Fin.ext_iff]

#print axioms target_means_observable
#print axioms targets_determine_law
#print axioms three_target_injective
#print axioms inversion_squared_bound
#print axioms activity_squared_bound
#print axioms binary_excess

end TimeFrequency.PredictiveActivity
