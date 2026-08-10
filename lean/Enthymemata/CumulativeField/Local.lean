import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Combinatorics.SimpleGraph.Basic
import Mathlib.Data.Fintype.Basic
import Mathlib.Data.Real.Basic
import Mathlib.Order.Monotone.Basic

/-!
# Cumulative fields — local calibration obligations

Candidate definitions for the recovered local construction. Their placement in
`Enthymemata` is deliberate: raw cumulative mass is only one possible local
calibration, and its density semantics remain under review.
-/

namespace CumulativeField

variable {V : Type*} [Fintype V]

/-- Cumulative incident coupling mass at vertex `u`, including the entire tie
class at `level`. -/
noncomputable def localCumulative (base : SimpleGraph V)
    (coupling : V → V → ℝ) (u : V) (level : ℝ) : ℝ := by
  classical
  exact (Finset.univ.filter
    (fun v => base.Adj u v ∧ coupling u v ≤ level)).sum
      (fun v => coupling u v)

/-- Directed edge score obtained by evaluating the source vertex's local
cumulative curve at the edge's own coupling. -/
noncomputable def localEdgeScore (base : SimpleGraph V)
    (coupling : V → V → ℝ) (u v : V) : ℝ :=
  localCumulative base coupling u (coupling u v)

/-- Pool all directed local cumulative curves at the same coupling level. This
is the exact algebraic aggregate against which an undirected global control can
later be compared after duplicate accounting is declared. -/
noncomputable def pooledCumulative (base : SimpleGraph V)
    (coupling : V → V → ℝ) (level : ℝ) : ℝ := by
  classical
  exact Finset.univ.sum (fun u => localCumulative base coupling u level)

/-- Whole-tie filtering makes a local edge score invariant under exchanging
two incident edges of equal coupling. -/
theorem localEdgeScore_eq_of_coupling_eq (base : SimpleGraph V)
    (coupling : V → V → ℝ) (u v w : V)
    (h : coupling u v = coupling u w) :
    localEdgeScore base coupling u v = localEdgeScore base coupling u w := by
  simp [localEdgeScore, h]

/-- Nonnegative couplings make every local cumulative curve monotone in its
level. This is the first open proof obligation for the candidate calibration. -/
theorem localCumulative_monotone (base : SimpleGraph V)
    (coupling : V → V → ℝ)
    (hnonneg : ∀ u v, base.Adj u v → 0 ≤ coupling u v) (u : V) :
    Monotone (localCumulative base coupling u) := by
  sorry

end CumulativeField
